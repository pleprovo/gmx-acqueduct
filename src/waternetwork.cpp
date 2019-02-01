#define _USE_MATH_DEFINES

#include "waternetwork.hpp"
#include "gromacs/pbcutil/pbc.h"

#include "boost/graph/graphviz.hpp"
#include "boost/graph/depth_first_search.hpp"
#include "boost/graph/bellman_ford_shortest_paths.hpp"

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
void add_bidirectional_edge(int u, int v, S s, Graph &g, bool direction = false)
{
    auto e1 = boost::add_edge(v, u, s, g).first;
    if (direction) {
	s.energy = 0.0;
    }
    auto e2 = boost::add_edge(u, v, s, g).first;
    g[e1].reverse_edge = e2;
    g[e2].reverse_edge = e1;
}


// void print_predecessor_path(std::vector<Graph::vertex_descriptor> parent,
// 			    Traits::vertex_descriptor u,
// 			    Traits::vertex_descriptor v)
// {
//     using path_t = std::vector<Graph::edge_descriptor>;
//     path_t path;    
//     for(Graph::vertex_descriptor u = g[v].predecessor; u != v; v=u, u=g[v].predecessor) {
//     	std::pair<Graph::edge_descriptor, bool> edge_pair = boost::edge(u,v,g);
//     	path.push_back( edge_pair.first );
//     }
        
//     std::cout << "Shortest Path from v1 to v6:" << std::endl;
//     for(path_t::reverse_iterator riter = path.rbegin(); riter != path.rend(); ++riter) {
//         Graph::vertex_descriptor u_tmp = boost::source(*riter, g);
//         Graph::vertex_descriptor v_tmp = boost::target(*riter, g);
//         Graph::edge_descriptor e_tmp = boost::edge(u_tmp, v_tmp, g).first;
	
//     	std::cout << "  " << g[u_tmp].id << " -> " << g[v_tmp].id << "    (weight: " << g[e_tmp].length << ")" << std::endl;
//     }
// }

void printPath(std::vector<Graph::vertex_descriptor> parent, int i, int j)
{
    // Base Case : If j is source
    if ( parent.at(j) == i )
        return;

    printPath(parent, i, parent.at(j));

    std::cout << parent.at(j) << " ";
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

double computeEnergy1(const double r, const double cosine,
		     const double r_on = 5.5,
		     const double r_off = 6.5,
		     const double theta_on = 0.25,
		     const double theta_off = 0.0301)
{
    static float C = 3855; /* epsilon*sigma^6*sqrt(2/3) */
    static float D = 738; /* epsilon*sigma^4*sqrt(2/3) */
    // static float theta_on = 0.25;
    // static float theta_off = 0.0301;
    // static float r_on = 5.5;
    // static float r_off = 6.5;
    double r_switch = switch_function(r, r_on, r_off);
    double theta_switch = 1-switch_function(pow(cosine, 2), theta_off, theta_on); 
    double energy = -((C/pow(r, 6.0))-(D/pow(r, 4.0)));
    energy *= pow(cosine,4.0);
    energy *= r_switch;
    energy *= theta_switch;
    return energy;
}

template<class EnergyMap, class DistanceMap, class AngleMap>
class edge_writer {
public:
    edge_writer(EnergyMap e, DistanceMap d, AngleMap a) : em(e),dm(d),am(a) {}
    template <class Edge>
    void operator()(std::ostream &out, const Edge& e) const {
	out << "[energy=\"" << em[e]
	    << "\", distance=\"" << dm[e]
	    << "\", angle=\"" << am[e]<< "\"]";
    }
private:
    EnergyMap em;
    DistanceMap dm;
    AngleMap am;
};

template<class EnergyMap, class DistanceMap, class AngleMap>
inline edge_writer<EnergyMap,DistanceMap,AngleMap> 
make_edge_writer(EnergyMap e,DistanceMap d,AngleMap a) {
    return edge_writer<EnergyMap,DistanceMap,AngleMap>(e,d,a);
}

WaterNetwork::WaterNetwork()
    : TrajectoryAnalysisModule("waternetwork", "Water network analysis tool"),
      alphavalue_(0.65), lengthon_(5.5), lengthoff_(6.5), angleon_(0.25), angleoff_(0.0301)
{
    WaterSearchPtr_ = std::make_shared<AlphaShapeSearch>();
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
    options->addOption(gmx::SelectionOption("alpha").store(&calpha_)
		       .description(""));    
    options->addOption(gmx::SelectionOption("protein").store(&protein_)
		       .description(""));
    options->addOption(gmx::SelectionOption("source").store(&source_)
		       .description("Define a group as a source for maximum flow analysis, must be used with -sink option"));
    options->addOption(gmx::SelectionOption("sink").store(&sink_)
		       .description("Define a group as a sink for maximum flow analysis, must be used with -source option"));


    options->addOption(gmx::DoubleOption("alphavalue").store(&alphavalue_)
		       .description("Alpha value for alpha shape computation (default 1.0nm)"));
    options->addOption(gmx::DoubleOption("don").store(&lengthon_)
		       .description("Distance Cutoff for energy calculation (default 0.25 nm)"));
    options->addOption(gmx::DoubleOption("doff").store(&lengthoff_)
		       .description("Distance Cutoff for energy calculation (default 0.25 nm)"));
    options->addOption(gmx::DoubleOption("aon").store(&angleon_)
		       .description("Angle cutoff for energy calculation (default cos(90))"));
    options->addOption(gmx::DoubleOption("aoff").store(&angleoff_)
		       .description("Angle cutoff for energy calculation (default cos(90))")); 

    
    //TODO link threshold energy to angle and distance
    //TODO write surface output
    //TODO write graph output
    
    settings->setFlag(gmx::TrajectoryAnalysisSettings::efRequireTop);
    settings->setFlag(gmx::TrajectoryAnalysisSettings::efUseTopX);
}

struct groupInfo {
    bool isDonor;
    int atom;
    int numberHydrogen;
};

void makeDonorAcceptorLists(gmx::Selection &selection, t_topology *top,
			    std::vector<std::pair<int, std::string> > &protIndices)
{
    std::vector<std::string> atomname{"NE", "NH1", "NH2", "ND1", "ND2", "OD1", "OD2", "NE1",
     	    "NE2", "OE1", "OE2", "NZ", "OG", "OG1", "OH", "O", "N", "SG"};
    unsigned int numAtoms = selection.posCount();
    gmx::ConstArrayRef<int> indices = selection.atomIndices();

    for (auto i : indices) {
    	std::string name(*top->atoms.atomname[i]);
        auto foo = std::find(names.begin(), names.end(), name);
	if (foo != names.end()) {
	    std::string plop(*top->atoms.atomname[i+1]);
	    protIndices.push_back(std::pair<int, std::string>(i, name));
	}
    }

    for (auto plop : protIndices) {
	if (plop.second == "SG") {
	    std::cout << plop.first << " "
		      << plop.second << " "
		      << *top->atoms.atomname[plop.first+1] <<  std::endl;
	}
	if (plop.second == "OG") {
	    std::cout << plop.first << " "
		      << plop.second << " "
		      << *top->atoms.atomname[plop.first+1] <<  std::endl;
	}
	if (plop.second == "N") {
	    std::cout << plop.first << " "
		      << plop.second << " "
		      << *top->atoms.atomname[plop.first+1] <<  std::endl;
	}
    }

    // loop over atoms and find potential donor or acceptor as in the atoms names vector
    // check is they are protonated by check if they are fullowed by hydrogen or correct number of hydrogen
    // write a struct containing all those guys and store it into a vector (include water ?) 
}


void WaterNetwork::initAnalysis(const gmx::TrajectoryAnalysisSettings &settings,
				const gmx::TopologyInformation        &top)
{
    /* Init Selection */

    solvent_.initOriginalIdsToGroup(top.topology(), INDEX_RES);
    nb_water_ = solvent_.posCount()/3;

    WaterSearchPtr_->setAlpha(alphavalue_);
    
    if (protein_.isValid()) {
	makeDonorAcceptorLists(protein_, top.topology(), protIndices_);
    }
    
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
    data_.setColumnCount(0, 5);

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
    outfilebinary = std::ofstream("file.binary", std::ios::out | std::ios::binary);
    outfilenormal = std::ofstream("file.dat", std::ios::out);

    /* Test the hydrogen bonds calculation */	
    // double distance = 3.0;
    // CGAL::Vector_3<K> OO{1.0, 0.0, 0.0};
    // CGAL::Vector_3<K> OH;
    // std::vector<double> results;
    // for (int i = 0; i < 100; i++) {
    // 	double angle = i*M_PI/200;
    // 	OH = CGAL::Vector_3<K>(cos(angle), sin(angle), 0.0);
    //     float cos_angle = OH*OO;
    // 	results.push_back(computeEnergy(distance, cos_angle));
    // 	std::cout << results.back() << " ";
    // }
    // std::cout << std::endl;
}


void
WaterNetwork::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
			   gmx::TrajectoryAnalysisModuleData *pdata)
{
    gmx::AnalysisDataHandle         dh     = pdata->dataHandle(data_);
    const gmx::Selection           &sourcesel = pdata->parallelSelection(source_);
    // const gmx::Selection           &proteinsel = pdata->parallelSelection(protein_);
    const gmx::Selection           &calphasel = pdata->parallelSelection(calpha_);
    const gmx::Selection           &sinksel = pdata->parallelSelection(sink_);
    const gmx::Selection           &watersel = pdata->parallelSelection(solvent_);

    /* Converstion of positions set to cgal point vectors */
    std::vector<Point_3> alphaPoints;
    // std::vector<Point_3> protVec = fromGmxtoCgalPosition<Point_3>(proteinsel.coordinates());
    std::vector<Point_3> watersVec = fromGmxtoCgalPosition<Point_3>(watersel.coordinates());
    std::vector<Point_3> oxygenVec = fromGmxtoCgalPosition<Point_3>(watersel.coordinates(), 3);
    std::vector<Point_3> sourceVec = fromGmxtoCgalPosition<Point_3>(sourcesel.coordinates());
    std::vector<Point_3> sinkVec = fromGmxtoCgalPosition<Point_3>(sinksel.coordinates());
    
    std::vector<int> buriedWaterVector;
    
    std::cout << " >>>> Frame Started " <<  frnr << std::endl;
    /* Alpha shape computation */
    /* Input : List of points*/
    /* Output : List of point in alpha shape or all points in initial selection*/
    if (calphasel.isValid()) {
	alphaPoints = fromGmxtoCgalPosition<Point_3>(calphasel.coordinates());
    	alphaShapeModulePtr_->build(alphaPoints);    
    	buriedWaterVector = alphaShapeModulePtr_->search(oxygenVec);
	// std::cout << " >> Alpha Shape done" << std::endl;
    } else {
	buriedWaterVector = std::vector<int>(nb_water_);
	std::iota(std::begin(buriedWaterVector), std::end(buriedWaterVector), 0);
    }
    
    /* Triangulation */
    /* Input : List of point */
    /* Output : Triangulation */   
    DelaunayWithInfo DT;
    Graph g;
    Graph gr;
    Graph gp;
    Ugraph gu;
    Graph subg = g;
    int count = 0;

    for (auto &indice : buriedWaterVector) {
	DelaunayWithInfo::Vertex_handle vh = DT.insert(oxygenVec.at(indice));
	vh->info() = Info{count, false, std::make_shared<Point_3>(watersVec.at(3*indice+1)),
			  std::make_shared<Point_3>(watersVec.at(3*indice+2))};
	boost::add_vertex(Atom {count}, g);
	boost::add_vertex(Atom {count}, gr);
	boost::add_vertex(Atom {count}, gu);
	boost::add_vertex(Atom {count}, gp);
	count++;
    }
    
    CGAL_assertion(DT.number_of_vertices() == oxygenVec.size());

    // for (auto &indice : protIndices_) {
    // 	DelaunayWithInfo::Vertex_handle vh = DT.insert(protVec.at(indice.first-protVec.size()));
    // 	vh->info() = Info{count, true};
    // 	boost::add_vertex(Atom {count}, g);
    // 	boost::add_vertex(Atom {count}, gr);
    // 	count++;	
    // }
    
    DelaunayWithInfo::Vertex_handle Source_handle = DT.insert(sourceVec.at(0));
    Source_handle->info() = Info{count, true};
    boost::add_vertex(Atom{count}, g);
    boost::add_vertex(Atom{count}, gr);
    boost::add_vertex(Atom {count}, gu);
    boost::add_vertex(Atom {count}, gp);
    count++;
    
    DelaunayWithInfo::Vertex_handle Sink_handle = DT.insert(sinkVec.at(0));
    Sink_handle->info() = Info{count, true};
    boost::add_vertex(Atom{count}, g);
    boost::add_vertex(Atom{count}, gr);
    boost::add_vertex(Atom {count}, gu);
    boost::add_vertex(Atom {count}, gp);

    // std::cout << " >> Triangulation done "<< DT.number_of_vertices()
    //  	      << " " << oxygenVec.size() << std::endl;
    
    /* Compute edge energy */
    /* Input : Triangulation */
    /* Output : Graph */
    CGAL::Vector_3<K> OO, OH11, OH12, OH21, OH22, OH1, OH2;
    double distance;
    int k, l;
    DelaunayWithInfo::Vertex_handle v1, v2;
    std::vector<std::pair<std::pair<int, int>, double> > plopVec;
    std::vector<double> hb_map;
    for(DelaunayWithInfo::Finite_edges_iterator ei=DT.finite_edges_begin();
	ei!=DT.finite_edges_end(); ++ei) {
	
	// std::cout << " ------------------ " << std::endl;
	v1 = ei->first->vertex(ei->second);
        v2 = ei->first->vertex(ei->third);
	// std::cout << v1->info().id << " " << v2->info().id << std::endl;
	// std::cout << v1->info().isSuperNode << " " << v2->info().isSuperNode << std::endl;

	HydrogenBond hb{distance, 0.0, 0.0};
	bool way = false;
	if (v1->info().isSuperNode or v2->info().isSuperNode) {
	    
	    // std::cout << "Super " << std::endl;
	    OO = v1->point()-v2->point();
	    distance = 10.0*CGAL::sqrt(OO*OO);
	    hb.angle = 1.0;
	    hb.energy = computeEnergy(distance, hb.angle);
	    
	} else {
	    std::vector<double> energies(4, 0.0);
	    // std::cout << "Normal " << std::endl;
	    k = buriedWaterVector.at(v1->info().id);
	    l = buriedWaterVector.at(v2->info().id);
	    OO = v1->point() - v2->point();
	    
	    distance = 10.0*CGAL::sqrt(OO*OO);

	    OH11 = *v1->info().h1 - v1->point();
	    OH12 = *v1->info().h2 - v1->point();
	    OH21 = *v2->info().h1 - v2->point();
	    OH22 = *v2->info().h2 - v2->point();
	    
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
		way = true;
	    }
	    hb.length = distance;
	    hb.angle = coss.at(ind);
	    hb.energy = *result;
	}

	if ( hb.energy != 0.0 ) {
	    // std::cout << "Distance : " << hb.length
	    // 	      << ". Angle : " << hb.angle
	    // 	      << ", Energy : " << hb.energy << std::endl;
	    if (way) {
		add_bidirectional_edge(v1->info().id, v2->info().id, hb, g, false);
		add_bidirectional_edge(v1->info().id, v2->info().id, hb, gr, true);
		plopVec.push_back(std::pair<std::pair<int, int>, double>(
				      std::pair<int, int>(v1->info().id, v2->info().id),
				      hb.energy));
		hb_map.push_back(hb.energy);
		hb.energy *= -1.0;
		add_bidirectional_edge(v1->info().id, v2->info().id, hb, gp, true);
	    } else {
		add_bidirectional_edge(v2->info().id, v1->info().id, hb, g, false);
		add_bidirectional_edge(v2->info().id, v1->info().id, hb, gr, true);
		plopVec.push_back(std::pair<std::pair<int, int>, double>(
				      std::pair<int, int>(v2->info().id, v1->info().id),
				      hb.energy));
		hb_map.push_back(hb.energy);
		hb.energy *= -1.0;
		add_bidirectional_edge(v1->info().id, v2->info().id, hb, gp, true);
	    }
	    
	    boost::add_edge(v1->info().id, v2->info().id, hb, gu);
	    std::cout << hb.energy << ", ";
	}
    }
    std::cout << std::endl;
    std::cout << " >> Energies done " << std::endl;
    
    /* Graph Analysis */
    double flow = do_max_flow(g, Source_handle->info().id, Sink_handle->info().id);    
    double flowf = do_max_flow(gr, Source_handle->info().id, Sink_handle->info().id);
    double flowb = do_max_flow(gr, Sink_handle->info().id, Source_handle->info().id);
    
    std::cout << " >> Flow done " << flow << " " << flowf << " " << flowb << std::endl;
    
    /* Connected Component */
    std::vector<int> component (boost::num_vertices (g));
    size_t num_components = boost::connected_components (gu, component.data()/*&component[0]*/);

    // if (component.at(Source_handle->info().id) != component.at(Sink_handle->info().id)) {
    // 	std::cout << " Number of Components : " << num_components 
    // 		  << ", Component source : " << component.at(Source_handle->info().id)
    // 		  << ", Component sink : " << component.at(Sink_handle->info().id) << std::endl;
    // }
    
    int component_size = 0;
    int component_ind = component.at(Source_handle->info().id);
    for (int i = 0; i < component.size(); i++) {
    	if (component.at(i) == component_ind) {
    	    // boost::add_vertex(Atom{i}, subg);
    	    // auto iter = boost::adjacent_vertices(i, g);
	    component_size++;
    	}
    }

    // for (auto &elem : dist) {
    // 	std::cout << elem << " "; 
    // }
    // std::cout << "\n";
    
    std::cout << " >> Components done " << component_size << std::endl;
    
    int N = num_vertices(gu);
    // Init pred vector 
    std::vector<std::size_t> pred(N);
    for (int i = 0; i< pred.size(); i++) {
	pred.at(i) = i;
    }

    // Init distance vector with source distance to zero
    std::vector<float> dist(N, (std::numeric_limits < short >::max)());
    dist[Source_handle->info().id] = 0;
    for (auto &elem : dist ) {
    	std::cout << elem << " ";
    }    
    std::cout << std::endl;
    
    bool r = boost::bellman_ford_shortest_paths
    	(gp, N, boost::weight_map(get(&HydrogenBond::energy, gp)).distance_map(&dist[0]).
    	 predecessor_map(&pred[0]).root_vertex(Sink_handle->info().id));
    if (r) {
	for (auto &elem : pred ) {
	    std::cout << elem << " ";
	}
	std::cout << std::endl;
	for (auto &elem : dist ) {
	    std::cout << elem << " ";
	}    
	std::cout << std::endl;
	std::cout << " >> Path from " << Source_handle->info().id << " to "
		  << Sink_handle->info().id << std::endl;
	printPath(pred, Source_handle->info().id, Sink_handle->info().id);
	std::cout << " >> Path done " << std::endl;
    } else {
	std::cout << " >> Path cycle problem " << std::endl;
    }
    /* Output Writing */
    
    // boost::write_graphviz(outfilebinary, gr,
    // 			  boost::make_label_writer(get(&Atom::id, gr)),
    // 			  make_edge_writer(get(&HydrogenBond::energy, gr),
    // 					   get(&HydrogenBond::length,gr),
    // 					   get(&HydrogenBond::angle, gr)));
    
    
    // boost::write_graphviz(outfilenormal, gr,
    // 			  boost::make_label_writer(get(&Atom::id, gr)),
    // 			  make_edge_writer(get(&HydrogenBond::energy, gr),
    // 					   get(&HydrogenBond::length,gr),
    // 					   get(&HydrogenBond::angle, gr)));
    
	
    // if (frnr == 100) {
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
    // 	if (calphasel.isValid()) {
    // 	    std::ofstream oss1;
    // 	    std::string ofs;
    // 	    oss1.open("surface.off");
    // 	    alphaShapeModulePtr_->writeOff(ofs);
    // 	    oss1 << ofs;
    // 	    oss1.close();
    // 	}
    // }

    
    
    dh.startFrame(frnr, fr.time);
    dh.setPoint(0, flow);
    dh.setPoint(1, flowf);
    dh.setPoint(2, flowb);
    dh.setPoint(3, 0.0 /*alphaShapeModulePtr_->volume()*/);
    dh.setPoint(3, component_size);
    dh.finishFrame();
    std::cout << " >> Output done" << std::endl; 
}


void WaterNetwork::finishAnalysis(int /*nframes*/)
{

}


void WaterNetwork::writeOutput()
{
    outfilebinary.close();
    outfilenormal.close();
}


int main(int argc, char *argv[])
{
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<WaterNetwork>(argc, argv);
}

