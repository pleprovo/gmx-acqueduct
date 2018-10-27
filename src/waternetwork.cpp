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

    if ( r_off > r && r > r_on ) {
	sw = pow(pow(r,2)-pow(r_off,2),2);
	sw *= pow(r_off,2)+2*pow(r,2)-3*pow(r_on,2);
	sw /= pow(pow(r_off,2)-pow(r_on,2),3);	
    } else if ( r <= r_on ) {
	sw = 1.0;
    }
    return sw;
}


double computeEnergy(const double r, const double cosine)
{
    static float C = 3855; /* epsilon*sigma^6*sqrt(2/3) */
    static float D = 738; /* epsilon*sigma^4*sqrt(2/3) */
    static float theta_on = 0.25;
    static float theta_off = 0.0301;
    static float r_on = 5.5;
    static float r_off = 6.5;
    double r_switch = switch_function(r, r_on, r_off);
    double theta_switch = 1-switch_function(pow(cosine, 2), theta_off, theta_on); 
    double energy = -((C/pow(r, 6.0))-(D/pow(r, 4.0)));
    energy *= pow(cosine,4.0);
    energy *= r_switch;
    energy *= theta_switch;
    return energy;
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
    options->addOption(gmx::SelectionOption("aplha").store(&calpha_)
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

    //makeDonorAcceptorLists(calpha_, top.topology());
    
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
    std::vector<Point_3> alphaPoints;
    std::vector<Point_3> watersVec = fromGmxtoCgalPosition<Point_3>(watersel.coordinates());
    std::vector<Point_3> oxygenVec = fromGmxtoCgalPosition<Point_3>(watersel.coordinates(), 3);
    std::vector<Point_3> sourceVec = fromGmxtoCgalPosition<Point_3>(sourcesel.coordinates());
    std::vector<Point_3> sinkVec = fromGmxtoCgalPosition<Point_3>(sinksel.coordinates());
    
    std::vector<int> buriedWaterVector;

    // std::cout << " ------------------ " << std::endl;
    // std::cout << " ------------------ " << std::endl;
    /* Alpha shape computation */
    /* Input : List of points*/
    /* Output : List of point in alpha shape or all points in initial selection*/
    if (false) {
	alphaPoints = fromGmxtoCgalPosition<Point_3>(calphasel.coordinates());
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

    for (auto &indice : buriedWaterVector) {
	DelaunayWithInfo::Vertex_handle vh = DT.insert(oxygenVec.at(indice));
	vh->info() = Info{count, false};
	boost::add_vertex(Atom {count}, g);
	count++;
    }
    
    CGAL_assertion(DT.number_of_vertices() == oxygenVec.size());
    
    // std::cout << DT.number_of_vertices() << " " << oxygenVec.size() << std::endl;

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
    CGAL::Vector_3<K> OO, OH11, OH12, OH21, OH22;
    double distance;
    int k, l;
    DelaunayWithInfo::Vertex_handle v1, v2;
    
    std::vector<std::pair<std::pair<int, int>, double> > plopVec;
    for(DelaunayWithInfo::Finite_edges_iterator ei=DT.finite_edges_begin();
	ei!=DT.finite_edges_end(); ++ei) {
	
	// std::cout << " ------------------ " << std::endl;
	v1 = ei->first->vertex(ei->second);
        v2 = ei->first->vertex(ei->third);
	// std::cout << v1->info().id << " " << v2->info().id << std::endl;
	// std::cout << v1->info().isSuperNode << " " << v2->info().isSuperNode << std::endl;

	HydrogenBond hb{distance, 0.0, 0.0};
	
	if (v1->info().isSuperNode or v2->info().isSuperNode) {
	    
	    // std::cout << "Super " << std::endl;
	    OO = v1->point()-v2->point();
	    distance = 10.0*CGAL::sqrt(OO*OO);
	    hb.angle = 1.0;
	    hb.energy = computeEnergy(distance, 1.0);
	    
	} else {
	    std::vector<double> energies(4, 0.0);
	    // std::cout << "Normal " << std::endl;
	    k = buriedWaterVector.at(v1->info().id);
	    l = buriedWaterVector.at(v2->info().id);
	    OO = oxygenVec.at(k) - oxygenVec.at(l);

	    distance = 10.0*CGAL::sqrt(OO*OO);
	
	    OH11 = watersVec.at(3*k+1) - oxygenVec.at(k);
	    OH12 = watersVec.at(3*k+2) - oxygenVec.at(k);
	    OH21 = watersVec.at(3*l+1) - oxygenVec.at(l);
	    OH22 = watersVec.at(3*l+2) - oxygenVec.at(l);

	    OO = OO / CGAL::sqrt(OO*OO);
	    OH11 = OH11 / CGAL::sqrt(OH11*OH11);
	    OH12 = OH12 / CGAL::sqrt(OH12*OH12);
	    OH21 = OH21 / CGAL::sqrt(OH21*OH21);
	    OH22 = OH22 / CGAL::sqrt(OH22*OH22);

	    std::vector<double> coss;
	    coss.push_back(OH11 * OO);
	    coss.push_back(OH12 * OO);
	    coss.push_back(OH21 * -OO);
	    coss.push_back(OH22 * -OO);
	    
	    // std::cout << "cos11 : " << coss.at(0) << ", "
	    // 	      << "cos12 : " << coss.at(1) << ", "
	    // 	      << "cos21 : " << coss.at(2) << ", "
	    // 	      << "cos22 : " << coss.at(3) << std::endl;
	    
	    for (unsigned int i = 0; i < coss.size(); i++) {
		if (coss.at(i) > 0.0) {
		    energies.at(i) = computeEnergy(distance, coss.at(i));
		    // std::cout << "energy" << i << " : " << energies.at(i) << std::endl;
		}
	    }
	    auto result = std::max_element(energies.begin(), energies.end());

	    int ind = std::distance(energies.begin(), result);
	    if ( ind < 2 ) {
		int temp;
	        temp = k;
		k = l;
		l = temp;		
	    }
	    hb.length = distance;
	    hb.angle = coss.at(ind);
	    hb.energy = *result;
	}

	if ( hb.energy > 0.001 ) {
	    // std::cout << "Distance : " << hb.length
	    // 	      << ". Angle : " << hb.angle
	    // 	      << ", Energy : " << hb.energy << std::endl;
	
	    add_bidirectional_edge(v1->info().id, v2->info().id, hb, g);
	    plopVec.push_back(std::pair<std::pair<int, int>, double>(
				  std::pair<int, int>(v1->info().id, v2->info().id),
				  hb.energy));
	}
    }
    // std::cout << " 0000000000000000000 " << std::endl;
    
    /* Graph Analysis */
    double flow = do_max_flow(g, Source_handle->info().id, Sink_handle->info().id);
    double flowr = do_max_flow(g, Sink_handle->info().id, Source_handle->info().id);
    double flows = flow - flowr;

    // std::cout << " 1111111111111111111 " << std::endl;
    
    /* Output Writing */
    if (frnr == 10) {
    	std::ofstream oss;
    	oss.open("triangulation.dat");
    	oss << boost::num_vertices(g) << " " << boost::num_edges(g) << "\n";
        for(DelaunayWithInfo::Vertex_iterator vi=DT.vertices_begin();
    	    vi!=DT.vertices_end(); vi++) {
    	    oss << vi->point() << "\n";
    	}
	
    	for (auto &elem : plopVec) {
    	    oss << elem.first.first << " " << elem.first.second << " " << elem.second << "\n";
    	}
    	oss << Source_handle->point() << "\n";
    	oss << Sink_handle->point() << "\n";
    	oss.close();

        // std::ofstream oss1;
    	// std::string ofs;
    	// oss1.open("surface.off");
    	// alphaShapeModulePtr_->writeOff(ofs);
    	// oss1 << ofs;
    	// oss1.close();
    }

    dh.startFrame(frnr, fr.time);
    dh.setPoint(0, flow);
    dh.setPoint(1, flowr);
    dh.setPoint(2, flows);
    dh.setPoint(3, 0.0 /*alphaShapeModulePtr_->volume()*/);
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

