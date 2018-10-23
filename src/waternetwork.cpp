#define _USE_MATH_DEFINES

#include "waternetwork.hpp"
#include "gromacs/pbcutil/pbc.h"

#include <iostream>
#include <fstream>
#include <limits>
#include <math.h>

struct custom_comparator {
    bool operator()(const std::pair<int, int>& a,
                    const std::pair<int, int>& b) const {
        return less_comparator(std::minmax(a.first, a.second),
                               std::minmax(b.first, b.second));
    }

    std::less<std::pair<int, int>> less_comparator;
};

/* Graph Stuff */
template <typename S>
void add_bidirectional_edge(int u, int v, S s, Graph &g)
{
    auto e1 = boost::add_edge(v, u, s, g).first;
    s.energy = 0.0;
    auto e2 = boost::add_edge(u, v, s, g).first;
    g[e1].reverse_edge = e2;
    g[e2].reverse_edge = e1;
}

void print_predecessor_path(Graph &g, Traits::vertex_descriptor v)
{
    using path_t = std::vector<Graph::edge_descriptor>;
    path_t path;    
    for(Graph::vertex_descriptor u = g[v].predecessor; u != v; v=u, u=g[v].predecessor) {
    	std::pair<Graph::edge_descriptor, bool> edge_pair = boost::edge(u,v,g);
    	path.push_back( edge_pair.first );
    }
        
    std::cout << "Shortest Path from v1 to v6:" << std::endl;
    for(path_t::reverse_iterator riter = path.rbegin(); riter != path.rend(); ++riter) {
        Graph::vertex_descriptor u_tmp = boost::source(*riter, g);
        Graph::vertex_descriptor v_tmp = boost::target(*riter, g);
        Graph::edge_descriptor e_tmp = boost::edge(u_tmp, v_tmp, g).first;
	
    	std::cout << "  " << g[u_tmp].id << " -> " << g[v_tmp].id << "    (weight: " << g[e_tmp].length << ")" << std::endl;
    }
}

double do_max_flow(Graph &g, const int source, const int sink)
{
    auto idx = get(&Atom::id, g);
    auto cap    = get(&HydrogenBond::energy, g);
    auto rescap = get(&HydrogenBond::residual_energy, g);
    auto rev = get(&HydrogenBond::reverse_edge, g);

    double flow = boost::boykov_kolmogorov_max_flow(g, cap, rescap, rev, idx, source, sink);
    return flow;
}

/* Gromacs to CGAL convertion function */    
template <class T>
std::vector<T> fromGmxtoCgalPosition(const gmx::ConstArrayRef<rvec> &coordinates,
				     const int increment=1)
{
    std::vector<T> cgalPositionVector;   
    for (unsigned int i = 0; i < coordinates.size(); i += increment) {
	cgalPositionVector.push_back(T(coordinates.at(i)[XX],
				       coordinates.at(i)[YY],
				       coordinates.at(i)[ZZ]));
    }
    return cgalPositionVector;
}


/* Hydrogen Bond stuff */
double switch_function(double r, double r_on, double r_off)
{
    double sw = 0.0;

    if ( r_off >= r && r >= r_on ) {
	sw = (pow(pow(r,2)-pow(r_off,2),2)*(pow(r_off,2)+2*pow(r,2)-3*pow(r_on,2)))/pow(pow(r_off,2)-pow(r_on,2),3);
    } else if ( r < r_on ) {
	sw = 1.0;
    }
    return sw;
}


double computeEnergy(const double r, const double cosine)
{
    static float C = 3855; /* epsilon*sigma^6*sqrt(2/3) */
    static float D = 738; /* epsilon*sigma^4*sqrt(2/3) */
    static float theta_on = 0.0669873;
    static float theta_off = 0.0;
    static float r_on = 5.5;
    static float r_off = 6.5;
    double energy = 0.0;
    double angle = std::acos(cosine);
    double r_switch = switch_function(r, r_on, r_off);
    double theta_switch = 1-switch_function(pow(cosine, 2), theta_off, theta_on); 
    energy = -((C/pow(r, 6.0))-(D/pow(r, 4.0)))*r_switch*theta_switch*pow(cosine, 4.0);
    return energy;
}


HydrogenBond computeHB(const rvec acceptor, const rvec donor, const rvec hydrogen)
{
    double energy = 0.0;
    rvec AD, DH;
    HydrogenBond out;
    rvec_sub(hydrogen, donor, DH);
    rvec_sub(acceptor, donor, AD);
    
    double r = 10.0*norm(AD);
    double angle = cos_angle(DH, AD);
    if ( cos_angle(DH, AD) < -0.52 ) {
	energy = computeEnergy(r, cos(angle));
    } else {
	energy = 0.0;
    }
    
    out.length = r;
    out.angle = angle;
    out.energy = energy;
    
    return out;
}


void checkHB(const gmx::AnalysisNeighborhoodPair &pair,
	     const gmx::Selection &waterSelection,
	     const gmx::ConstArrayRef<int> &mappedIds,
	     std::vector<HydrogenBond> &edgeVector)
{
    auto donor = waterSelection.position(mappedIds.at(pair.testIndex())).x();
    auto hydrogen1 = waterSelection.position(mappedIds.at(pair.testIndex())+1).x();
    auto hydrogen2 = waterSelection.position(mappedIds.at(pair.testIndex())+2).x();
    auto acceptor = waterSelection.position(mappedIds.at(pair.refIndex())).x();

    HydrogenBond out1, out2, out;
    out1 = computeHB(acceptor, donor, hydrogen1);
    out2 = computeHB(acceptor, donor, hydrogen2);
    if (out1.energy > out2.energy) {
	out = out1;
    } else {
	out = out2;
    }
    
    
    edgeVector.push_back(out);							  
}


WaterNetwork::WaterNetwork()
    : TrajectoryAnalysisModule("waternetwork", "Water network analysis tool"),
      cutoff_(0.65)
{
    alphaShapeModulePtr_ = std::make_shared<AlphaShapeModule>();
    registerAnalysisDataset(&data_, "avedist");
}


void WaterNetwork::initOptions(gmx::Options                    *options,
			       gmx::TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "This is a template for writing your own analysis tools for",
        "GROMACS. The advantage of using GROMACS for this is that you",
        "have access to all information in the topology, and your",
        "program will be able to handle all types of coordinates and",
        "trajectory files supported by GROMACS. In addition,",
        "you get a lot of functionality for free from the trajectory",
        "analysis library, including support for flexible dynamic",
        "selections. Go ahead an try it![PAR]",
        "To get started with implementing your own analysis program,",
        "follow the instructions in the README file provided.",
        "This template implements a simple analysis programs that calculates",
        "average distances from a reference group to one or more",
        "analysis groups."
    };

    options->setDescription(desc);

    options->addOption(gmx::FileNameOption("o").filetype(gmx::eftPlot).outputFile()
		       .store(&fnDist_).defaultBasename("avedist")
		       .description("Collection of analysis properties through time"));
    options->addOption(gmx::SelectionOption("select").store(&solvent_).required()
		       .defaultSelectionText("Water")
		       .description("Groups to calculate graph properties (default Water)"));    
    options->addOption(gmx::SelectionOption("aplha").store(&calpha_).required()
		       .description(""));
    options->addOption(gmx::SelectionOption("source").store(&source_)
		       .description("Define a group as a source for maximum flow analysis, must be used with -sink option"));
    options->addOption(gmx::SelectionOption("sink").store(&sink_)
		       .description("Define a group as a sink for maximum flow analysis, must be used with -source option"));    
    options->addOption(gmx::DoubleOption("cutoff").store(&cutoff_)
		       .description("Cutoff for distance calculation (default 0.25 nm)"));

    settings->setFlag(gmx::TrajectoryAnalysisSettings::efRequireTop);
    settings->setFlag(gmx::TrajectoryAnalysisSettings::efUseTopX);
}


void makeDonorAcceptorLists(gmx::Selection &selection, t_topology *top)
{    
    // std::vector<int> donorIndices;
    // std::vector<int> acceptorIndices;
    // std::cout << "Number of atoms : " << numAtoms << std::endl;
    // std::cout << "Number of residues : " << numGroup << std::endl;
    // for (unsigned int i = 0; i < numAtoms; i++) {
    // 	std::string name(*top->atoms.atomname[i]);
    // 	if (name == "O") {
    // 	    acceptorIndices.push_back(i);
    // 	}
    // 	if (name == "N") {
    // 	    donorIndices.push_back(i);
    // 	}
    // 	if (name == "SG") {
    // 	    donorIndices.push_back(i);
    // 	}
    // 	if (name == "NZ") {
    // 	    donorIndices.push_back(i);
    // 	}
    // }
    // std::cout << "Number of donor :" << donorIndices.size() << std::endl;
    // std::cout << "Number of acceptor :" << acceptorIndices.size() << std::endl;
    
    // std::cout << std::endl;
}


void WaterNetwork::initAnalysis(const gmx::TrajectoryAnalysisSettings &settings,
				const gmx::TopologyInformation        &top)
{
    /* Init Selection */

    solvent_.initOriginalIdsToGroup(top.topology(), INDEX_RES);
    nb_water_ = solvent_.posCount()/3;

    makeDonorAcceptorLists(calpha_, top.topology());
    
    /* Init neihborsearch cutoff value */
    // nb1_.setCutoff(cutoff_);
    // nb2_.setCutoff(cutoff_);
    // nb1_.setMode(gmx::AnalysisNeighborhood::SearchMode::eSearchMode_Grid);
    // nb2_.setMode(gmx::AnalysisNeighborhood::SearchMode::eSearchMode_Grid);
    
    /* Create ConstArray of indices for oxygen atoms */
    oxygenIndices_.reserve(nb_water_);
    for (int i = 0; i < nb_water_; i++) {
    	oxygenIndices_.push_back(3*i);
    }

    /* Test flow directionallity */

    
    
    
    /* Set the number of column to store time dependent data */
    data_.setColumnCount(0, 4);

    /* Init the average module  */ 
    avem_.reset(new gmx::AnalysisDataAverageModule());
    data_.addModule(avem_);

    /* Init the Plot module for the time dependent data */
    if (!fnDist_.empty()) {
	gmx::AnalysisDataPlotModulePointer plotm(
	    new gmx::AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnDist_);
        plotm->setTitle("Graph Statistics");
        plotm->setXAxisIsTime();
        //plotm->setYLabel("Distance (nm)");
        data_.addModule(plotm);
    }
}


void
WaterNetwork::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
			   gmx::TrajectoryAnalysisModuleData *pdata)
{
    gmx::AnalysisDataHandle         dh     = pdata->dataHandle(data_);
    const gmx::Selection           &sourcesel = pdata->parallelSelection(source_);
    const gmx::Selection           &calphasel = pdata->parallelSelection(calpha_);
    const gmx::Selection           &sinksel = pdata->parallelSelection(sink_);
    const gmx::Selection           &watersel = pdata->parallelSelection(solvent_);

    /* Converstion of positions set to cgal point vectors */
    std::vector<Point_3> alphaPoints = fromGmxtoCgalPosition<Point_3>(calphasel.coordinates());
    std::vector<Point_3> watersVec = fromGmxtoCgalPosition<Point_3>(watersel.coordinates());
    std::vector<Point_3> oxygenVec = fromGmxtoCgalPosition<Point_3>(watersel.coordinates(), 3);
    std::vector<Point_3> sourceVec = fromGmxtoCgalPosition<Point_3>(sourcesel.coordinates());
    std::vector<Point_3> sinkVec = fromGmxtoCgalPosition<Point_3>(sinksel.coordinates());
    
    std::vector<int> buriedWaterVector;
    
    /* Alpha shape computation */
    /* Input : List of points*/
    /* Output : List of point in alpha shape or all points in initial selection*/
    if (true) {
	alphaShapeModulePtr_->build(alphaPoints, 1.0);    
	buriedWaterVector = alphaShapeModulePtr_->locate(oxygenVec, true);
    } else {
	buriedWaterVector = std::vector<int>(nb_water_);
	std::iota(std::begin(buriedWaterVector), std::end(buriedWaterVector), 0);
    }

    /* Triangulation */
    /* Input : List of point */
    /* Output : Triangulation */   
    DelaunayWithInfo DT;
    Graph g;
    int count = 0;
    HydrogenBond hb;

    for (auto &indice : buriedWaterVector)
    {
	DelaunayWithInfo::Vertex_handle vh = DT.insert(oxygenVec.at(indice));
	vh->info() = Info{count, false};
	boost::add_vertex(Atom{count}, g);
	count++;
    }
    
    CGAL_assertion(DT.number_of_vertices() == oxygenVec.size());
        
    DelaunayWithInfo::Vertex_handle Source_handle = DT.insert(sourceVec.at(0));
    Source_handle->info() = Info{count, true};
    boost::add_vertex(Atom{count}, g);
    count++;
    DelaunayWithInfo::Vertex_handle Sink_handle = DT.insert(sinkVec.at(0));
    Sink_handle->info() = Info{count, true};
    boost::add_vertex(Atom{count}, g);


    /* Compute edge enrgy */
    /* Input : Triangulation */
    /* Output : Graph */
    CGAL::Vector_3<K> OOp, OOd, OH11, OH12, OH21, OH22;
    double distance;
    int k, l;
    DelaunayWithInfo::Vertex_handle v1, v2;
    double cos11, cos12, cos21, cos22;
    std::vector<double> energies(4, 0.0);
    for(DelaunayWithInfo::Finite_edges_iterator ei=DT.finite_edges_begin();
	ei!=DT.finite_edges_end(); ++ei) {
	v1 = ei->first->vertex(ei->second);
        v2 = ei->first->vertex(ei->third);

	if (v1->info().isSuperNode or v2->info().isSuperNode) {
	    OOp = v1->point()-v2->point();
	    distance = 10.0*CGAL::sqrt(OOp*OOp);
	    hb.energy = computeEnergy(distance, 1.0);
	    hb.length = distance;
	} else {
	    k = buriedWaterVector.at(v1->info().id);
	    l = buriedWaterVector.at(v2->info().id);

	    OOp = oxygenVec.at(k) - oxygenVec.at(l);
	    OOd = oxygenVec.at(l) - oxygenVec.at(k);

	    distance = 10.0*CGAL::sqrt(OOp*OOp);
	
	    OH11 = watersVec.at(3*k+1) - oxygenVec.at(k);
	    OH12 = watersVec.at(3*k+2) - oxygenVec.at(k);
	    OH21 = watersVec.at(3*l+1) - oxygenVec.at(l);
	    OH22 = watersVec.at(3*l+2) - oxygenVec.at(l);

	    OOp = OOp / CGAL::sqrt(OOp*OOp);
	    OOd = OOd / CGAL::sqrt(OOd*OOd);
	    OH11 = OH11 / CGAL::sqrt(OH11*OH11);
	    OH12 = OH12 / CGAL::sqrt(OH12*OH12);
	    OH21 = OH21 / CGAL::sqrt(OH21*OH21);
	    OH22 = OH22 / CGAL::sqrt(OH22*OH22);
	
	    cos11 = OH11 * OOp;
	    cos12 = OH12 * OOp;
	    cos21 = OH21 * OOd;
	    cos22 = OH22 * OOd;
	
	    energies.at(0) = computeEnergy(distance, cos11);
	    energies.at(1) = computeEnergy(distance, cos12);
	    energies.at(2) = computeEnergy(distance, cos21);
	    energies.at(3) = computeEnergy(distance, cos22);

	    auto result = std::max_element(energies.begin(), energies.end());
	    
	    hb.energy = *result;
	    hb.length = distance;

	    int ind = std::distance(energies.begin(), result);
	    if ( ind < 2 ) {
		int temp;
	        temp = k;
		k = l;
		l = temp;		
	    }	       
	}
	
	if (hb.energy != 0.0) {	
	    add_bidirectional_edge(v1->info().id, v2->info().id, hb, g);
	}
    }

    /* Graph Analysis */
    double flow = do_max_flow(g, Source_handle->info().id, Sink_handle->info().id);
    double flowr = do_max_flow(g, Sink_handle->info().id, Source_handle->info().id);
    double flows = flow + flowr;
    int iflow=0, iflowr=0, iflows=0, iflowp=0;
    if (flow > 0.0) iflow = 1;
    if (flowr > 0.0) iflow = -1;
    if (flow > 0.0 and flowr > 0.0) iflow = 2;
    /* Output Writing */

    
    // if (frnr == 10) {
    // 	std::ofstream oss;
    // 	oss.open("triangulation.dat");
    // 	oss << boost::num_vertices(g) << " " << boost::num_edges(g) << "\n";
    //     for(DelaunayWithInfo::Vertex_iterator vi=DT.vertices_begin();
    // 	    vi!=DT.vertices_end(); vi++) {
    // 	    oss << vi->point() << "\n";
    // 	}
	
    // 	for (auto &elem : plopVec) {
    // 	    oss << elem.first.first << " " << elem.first.second << " " << elem.second << "\n";
    // 	}
    // 	oss << Source_handle->point() << "\n";
    // 	oss << Sink_handle->point() << "\n";
    // 	oss.close();

    //     std::ofstream oss1;
    // 	std::string ofs;
    // 	oss1.open("surface.off");
    // 	alphaShapeModulePtr_->writeOff(ofs);
    // 	oss1 << ofs;
    // 	oss1.close();
    // }
    
    // std::cout << "\nBuried : " << buriedWaterVector.size() << " "
    // 	      << "DT(" << DT.number_of_vertices() << ", " << DT.number_of_edges() << ") "
    // 	      << "Gr(" << boost::num_vertices(g) << ", " << boost::num_edges(g) << ", "
    // 	      << "flow : " << flow << ")\n";
    
    /* Find all hydrogen bonds */

    // gmx::ConstArrayRef<int> oxygenArrayRef(oxygenIndices_);

    // std::set<int> sourceEdges;
    // std::set<int> sinkEdges;
    // std::vector<HydrogenBond> HBVector;
    
    // gmx::AnalysisNeighborhoodPair pair;
    // gmx::AnalysisNeighborhoodPositions sourcePos(sourcesel);
    // gmx::AnalysisNeighborhoodPositions calphaPos(calphasel);
    // gmx::AnalysisNeighborhoodPositions sinkPos(sinksel);    
    // gmx::AnalysisNeighborhoodPositions waterPos(watersel);
    
    // gmx::AnalysisNeighborhoodPositions oxygenPos = waterPos.indexed(oxygenArrayRef);
    
    // gmx::AnalysisNeighborhoodSearch searchAlpha = nb1_.initSearch(pbc, oxygenPos);
    // gmx::AnalysisNeighborhoodPairSearch pairSearchAlpha = searchAlpha.startPairSearch(calphaPos);
    // std::set<int> neighborSet;
    // std::vector<int> neighborVec;
    // while (pairSearchAlpha.findNextPair(&pair))
    // {
    //     neighborSet.insert(oxygenArrayRef.at(pair.refIndex()));
    // }
    
    // std::copy(neighborSet.begin(), neighborSet.end(), std::back_inserter(neighborVec));
    // gmx::ConstArrayRef<int> neighborArrayRef(neighborVec);
  
    // gmx::AnalysisNeighborhoodPositions neighborPos = waterPos.indexed(neighborArrayRef);
    // gmx::AnalysisNeighborhoodSearch search1 = nb1_.initSearch(pbc, neighborPos);
    // gmx::AnalysisNeighborhoodPairSearch pairSearch = search1.startPairSearch(neighborPos);
	
    // gmx::AnalysisNeighborhoodSearch search2 = nb2_.initSearch(pbc, oxygenPos);    
    // gmx::AnalysisNeighborhoodPairSearch pairSearchSource = search2.startPairSearch(sourcePos);
    // gmx::AnalysisNeighborhoodPairSearch pairSearchSink = search2.startPairSearch(sinkPos);

    // gmx::AnalysisNeighborhoodSearch search1 = nb1_.initSearch(pbc, oxygenPos);
    // gmx::AnalysisNeighborhoodPairSearch pairSearch = search1.startPairSearch(oxygenPos);
     
    // Graph g(nb_water_+2);

    // std::vector<std::pair<std::pair<int,int>, double>> hbn;
    
    // while (pairSearch.findNextPair(&pair))
    // {
    // 	if (pair.refIndex() != pair.testIndex())
    // 	{
    // 	    checkHB(pair, watersel, oxygenArrayRef, HBVector);
    // 	}
    // }
    
    // while (pairSearchSource.findNextPair(&pair))
    // {
    // 	sourceEdges.insert(pair.refIndex());
    // }

    // while (pairSearchSink.findNextPair(&pair))
    // {
    // 	sinkEdges.insert(pair.refIndex());
    // }
  
    dh.startFrame(frnr, fr.time);
    dh.setPoint(0, iflow);
    dh.setPoint(1, 0.0);
    dh.setPoint(2, 0.0);
    dh.setPoint(3, alphaShapeModulePtr_->volume());
    dh.finishFrame();
}


void WaterNetwork::finishAnalysis(int /*nframes*/)
{

}


void WaterNetwork::writeOutput()
{
    
}


int main(int argc, char *argv[])
{
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<WaterNetwork>(argc, argv);
}

