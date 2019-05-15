#define _USE_MATH_DEFINES

#include "waternetwork.hpp"
#include "gromacs/pbcutil/pbc.h"

#include "AnalysisAlphaShape.hpp"
#include "AnalysisNeighbors.hpp"
// #include "graph_module.hpp"

// #include "boost/graph/graphviz.hpp"
// #include "boost/graph/depth_first_search.hpp"
// #include "boost/graph/bellman_ford_shortest_paths.hpp"

#include <iostream>
#include <fstream>
#include <limits>
#include <math.h>

/* Graph Stuff */

// template <typename S>
// void add_bidirectional_edge(int u, int v, S s, Graph &g, bool direction = false)
// {
//     auto e1 = boost::add_edge(v, u, s, g).first;
//     if (direction) {
// 	s.energy = 0.0;
//     }
//     auto e2 = boost::add_edge(u, v, s, g).first;
//     g[e1].reverse_edge = e2;
//     g[e2].reverse_edge = e1;
// }
// */

// // void print_predecessor_path(std::vector<Graph::vertex_descriptor> parent,
// // 			    Traits::vertex_descriptor u,
// // 			    Traits::vertex_descriptor v)
// // {
// //     using path_t = std::vector<Graph::edge_descriptor>;
// //     path_t path;    
// //     for(Graph::vertex_descriptor u = g[v].predecessor; u != v; v=u, u=g[v].predecessor) {
// //     	std::pair<Graph::edge_descriptor, bool> edge_pair = boost::edge(u,v,g);
// //     	path.push_back( edge_pair.first );
// //     }
        
// //     std::cout << "Shortest Path from v1 to v6:" << std::endl;
// //     for(path_t::reverse_iterator riter = path.rbegin(); riter != path.rend(); ++riter) {
// //         Graph::vertex_descriptor u_tmp = boost::source(*riter, g);
// //         Graph::vertex_descriptor v_tmp = boost::target(*riter, g);
// //         Graph::edge_descriptor e_tmp = boost::edge(u_tmp, v_tmp, g).first;
	
// //     	std::cout << "  " << g[u_tmp].id << " -> " << g[v_tmp].id << "    (weight: " << g[e_tmp].length << ")" << std::endl;
// //     }
// // }
//  /*
// void printPath(std::vector<Graph::vertex_descriptor> parent, int i, int j)
// {
//     // Base Case : If j is source
//     if ( parent.at(j) == i )
//         return;

//     printPath(parent, i, parent.at(j));

//     std::cout << parent.at(j) << " ";
// }

// double do_max_flow(Graph &g, const int source, const int sink)
// {
//     auto idx = get(&Atom::id, g);
//     auto cap    = get(&HydrogenBond::energy, g);
//     auto rescap = get(&HydrogenBond::residual_energy, g);
//     auto rev = get(&HydrogenBond::reverse_edge, g);

//     double flow = boost::boykov_kolmogorov_max_flow(g, cap, rescap, rev, idx, source, sink);
//     return flow;
// }


int makeSites(gmx::Selection &selection, t_topology *top, std::vector<Site> &sites)
{
    std::vector<std::string> atomNames{
	"NE",
	    "NH1",
	    "NH2",
	    "ND1",
	    "ND2",
	    "OD1",
	    "OD2",
	    "NE1",
     	    "NE2",
	    "OE1",
	    "OE2",
	    "NZ",
	    "OG",
	    "OG1",
	    "OH",
	    "O",
	    "N",
	    "SG",
	    "OW"};

    gmx::ConstArrayRef<int> indices = selection.atomIndices();
    gmx::ConstArrayRef<int> resIndices = selection.mappedIds();
    for (int i = 0; i < indices.size(); i++)
    {
	int index = i;
    	std::string name(*top->atoms.atomname[index]);
        auto foo = std::find(atomNames.begin(), atomNames.end(), name);

	if (foo != atomNames.end())
	{
	    Site site;
	    site.index = index;
	    site.resIndex = resIndices.at(i);
	    site.nbHydrogen = 0;
	    int j = i+1;
	    bool flag = true;
	    while (flag)
	    {
	    	if (j >= indices.size())
	    	{
	    	    sites.push_back(site);
	    	    return sites.size();
	    	}
	    	if (std::string(*top->atoms.atomname[j+site.nbHydrogen])[0] == 'H')
	    	{
	    	    site.nbHydrogen++;
	    	    j++;
	    	}
	    	else
	    	{
	    	    flag=false;
	    	}
	    }
	    sites.push_back(site);
	}
    }
    return sites.size();
}


WaterNetwork::WaterNetwork()
    : TrajectoryAnalysisModule("waternetwork", "Water network analysis tool"),
      alphavalue_(0.65), lengthon_(5.5), lengthoff_(6.5), angleon_(0.25), angleoff_(0.0301)
{
    strategy_ = std::shared_ptr<AnalysisInterface>(new StrategyAlpha);
    //strategy_ = std::shared_ptr<AnalysisInterface>(new AnalysisNeighbors);
    registerAnalysisDataset(&filterData_, "filter");
    registerAnalysisDataset(&graphData_, "graph");
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

    options->addOption(gmx::FileNameOption("of").filetype(gmx::eftPlot).outputFile()
		       .store(&fnFilter_).defaultBasename("filter")
		       .description("Collection of analysis properties through time"));

    options->addOption(gmx::FileNameOption("og").filetype(gmx::eftPlot).outputFile()
		       .store(&fnGraph_).defaultBasename("graph")
		       .description("Collection of analysis properties through time"));
    
    options->addOption(gmx::SelectionOption("select").store(&solvent_).required()
		       .defaultSelectionText("SOL")
		       .description("Groups to calculate graph properties (default Water)"));    
    
    options->addOption(gmx::SelectionOption("protein").store(&protein_).required()
		       .defaultSelectionText("Protein")
		       .description(""));
    options->addOption(gmx::SelectionOption("alpha").store(&calpha_).required()
		       .defaultSelectionText("Calpha")
		       .description(""));    
    options->addOption(gmx::SelectionOption("source").store(&source_)
		       .description("Define a group as a source for maximum flow analysis, must be used with -sink option"));
    options->addOption(gmx::SelectionOption("sink").store(&sink_)
		       .description("Define a group as a sink for maximum flow analysis, must be used with -source option"));

    const char * const  allowed[] = { "alphashape", "neighbor"};
    std::string  str;
    int          filterType;
    options->addOption(gmx::StringOption("type").enumValue(allowed).store(&str)
		       .defaultValue("alphashape")
		       .storeEnumIndex(&filterType));

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
    
    settings->setFlag(gmx::TrajectoryAnalysisSettings::efRequireTop);
    settings->setFlag(gmx::TrajectoryAnalysisSettings::efUseTopX);

    /* Display option */
    std::cout << "Use Filter type : " << str << std::endl;
}


void WaterNetwork::initAnalysis(const gmx::TrajectoryAnalysisSettings &settings,
				const gmx::TopologyInformation        &top)
{
    /* Make Solvent Sites */
    solvent_.initOriginalIdsToGroup(top.topology(), INDEX_RES);
    std::vector<Site> solventSites;
    if (solvent_.isValid()) {
    	if (makeSites(solvent_, top.topology(), solventSites_) > 0)
    	{
    	    std::clog << "Solvent has : " << solventSites_.size()
		      << " residues."<< std::endl;
    	}
    }
    
    /* Make Protein Sites */
    protein_.initOriginalIdsToGroup(top.topology(), INDEX_RES);
    std::vector<Site> proteinSites;
    if (protein_.isValid()) {
	if (makeSites(protein_, top.topology(), proteinSites_) > 0)
	{
	    std::clog << "Protein has : " << proteinSites_.size()
		      << " residues." << std::endl;
	}
    }
    
    /* Filter Type for Analysis */
    
    /* Set the number of column to store time dependent data */
    filterData_.setColumnCount(0, 2); 
    graphData_.setColumnCount(0, 3);

    /* Init the average module  */ 
    avem_.reset(new gmx::AnalysisDataAverageModule());
    filterData_.addModule(avem_);

    /* Init the Plot module for the time dependent data */
    /* Solvent filtering data */
    if (!fnFilter_.empty()) {
	gmx::AnalysisDataPlotModulePointer plotm(
	    new gmx::AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnFilter_);
        plotm->setTitle("Filter Statistics");
        plotm->setXAxisIsTime();
        //plotm->setYLabel("Distance (nm)");
        filterData_.addModule(plotm);
    }

    /* Init the Plot module for the time dependent data */
    /* Network Data */
    if (!fnGraph_.empty()) {
	gmx::AnalysisDataPlotModulePointer plotm(
	    new gmx::AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnGraph_);
        plotm->setTitle("Graph Statistics");
        plotm->setXAxisIsTime();
        //plotm->setYLabel("Distance (nm)");
        graphData_.addModule(plotm);
    }
    
}


void
WaterNetwork::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
			   gmx::TrajectoryAnalysisModuleData *pdata)
{
    gmx::AnalysisDataHandle         filterData = pdata->dataHandle(filterData_);
    gmx::AnalysisDataHandle         graphData = pdata->dataHandle(graphData_);
    const gmx::Selection           &proteinsel = pdata->parallelSelection(protein_);
    const gmx::Selection           &calphasel = pdata->parallelSelection(calpha_);
    const gmx::Selection           &solventsel = pdata->parallelSelection(solvent_);

    
    /*
     * Store frame info
     */
    
    Frame frame{proteinsel, proteinSites_,
	    solventsel, solventSites_, calphasel};

    /*
     * Analyse the frame
     */
    Results results = strategy_->execute(frame);
    
    /* Filter Solvent Sites */
    /* Choose Filter */
    // int filtered = filterSolventSite<FilterStrategy>(solventSites, filterSelection);
    
    /* Concatenate Protein and Solvent Sites */
    // sites.insert(sites.end(),
    // 		 std::make_move_iterator(proteinSites.begin()),
    // 		 std::make_move_iterator(proteinSites.end()));
    // sites.insert(sites.end(),
    // 		 std::make_move_iterator(waterSites.begin()),
    // 		 std::make_move_iterator(waterSites.end()));
    
    /* Make Edges from Sites */
    // std::vector<Edges> makeEdges<EdgeStrategy>(Sites);

    /* Choose Edge evaluator */
    // std::property_map<Property> evaluateHydrogenbond(Edges);
    
    /* Build Graph  */
    // int success = makeGraph(Sites, Edges, property_map);

    /* Analyse Graph */
    /* Choose type of analysis */

    /* Write Graph */
    /* Find a way to write it efficiently */
   
    // Graph g;
    // Graph gr;
    // Graph gp;
    // Ugraph gu;
    // Graph subg = g;
    // */
    // int count = 0;

    // // if (proteinsel.isValid()) {
    // // 	for ()
    // // 	DelaunayWithInfo::Vertex_handle vh = DT.insert(oxygenVec.at(indice));
    // // 	vh->info() = Info{count, false, std::make_shared<Point_3>(watersVec.at(3*indice+1)),
    // // 			  std::make_shared<Point_3>(watersVec.at(3*indice+2))};
    // // 	boost::add_vertex(Atom {count}, g);
    // // 	boost::add_vertex(Atom {count}, gr);
    // // 	boost::add_vertex(Atom {count}, gu);
    // // 	count++;
    // // }
    
    // for (auto &indice : buriedWaterVector) {
    // 	DelaunayWithInfo::Vertex_handle vh = DT.insert(oxygenVec.at(indice));
    // 	vh->info() = Info{count, false, std::make_shared<Point_3>(watersVec.at(3*indice+1)),
    // 			  std::make_shared<Point_3>(watersVec.at(3*indice+2))};
    // 	/*
    // 	boost::add_vertex(Atom {count}, g);
    // 	boost::add_vertex(Atom {count}, gr);
    // 	boost::add_vertex(Atom {count}, gu);
    // 	boost::add_vertex(Atom {count}, gp);
    // 	count++;
    // }
    
    // CGAL_assertion(DT.number_of_vertices() == oxygenVec.size());

    // // for (auto &indice : protIndices_) {
    // // 	DelaunayWithInfo::Vertex_handle vh = DT.insert(protVec.at(indice.first-protVec.size()));
    // // 	vh->info() = Info{count, true};
    // // 	boost::add_vertex(Atom {count}, g);
    // // 	boost::add_vertex(Atom {count}, gr);
    // // 	count++;	
    // // }
    
    // std::cout << " >> Triangulation done "<< DT.number_of_vertices()
    //  	      << " " << oxygenVec.size() << std::endl;
    
    // DelaunayWithInfo::Vertex_handle Source_handle = DT.insert(sourceVec.at(0));
    // Source_handle->info() = Info{count, true};
    // boost::add_vertex(Atom{count}, g);
    // boost::add_vertex(Atom{count}, gr);
    // boost::add_vertex(Atom {count}, gu);
    // boost::add_vertex(Atom {count}, gp);
    // count++;
    
    // DelaunayWithInfo::Vertex_handle Sink_handle = DT.insert(sinkVec.at(0));
    // Sink_handle->info() = Info{count, true};
    // boost::add_vertex(Atom{count}, g);
    // boost::add_vertex(Atom{count}, gr);
    // boost::add_vertex(Atom {count}, gu);
    // boost::add_vertex(Atom {count}, gp);

    // // std::cout << " >> Triangulation done "<< DT.number_of_vertices()
    // //  	      << " " << oxygenVec.size() << std::endl;
    
    // /* Compute edge energy */
    // /* Input : Triangulation */
    // /* Output : Graph */
    // CGAL::Vector_3<K> OO, OH11, OH12, OH21, OH22, OH1, OH2;
    // double distance;
    // int k, l;
    // DelaunayWithInfo::Vertex_handle v1, v2;
    // std::vector<std::pair<std::pair<int, int>, double> > plopVec;
    // std::vector<double> hb_map;
    // for(DelaunayWithInfo::Finite_edges_iterator ei=DT.finite_edges_begin();
    // 	ei!=DT.finite_edges_end(); ++ei) {
	
    // 	// std::cout << " ------------------ " << std::endl;
    // 	v1 = ei->first->vertex(ei->second);
    //     v2 = ei->first->vertex(ei->third);
    // 	// std::cout << v1->info().id << " " << v2->info().id << std::endl;
    // 	// std::cout << v1->info().isSuperNode << " " << v2->info().isSuperNode << std::endl;

    // 	HydrogenBond hb{distance, 0.0, 0.0};
    // 	bool way = false;
    // 	if (v1->info().isSuperNode or v2->info().isSuperNode) {
	    
    // 	    // std::cout << "Super " << std::endl;
    // 	    OO = v1->point()-v2->point();
    // 	    distance = 10.0*CGAL::sqrt(OO*OO);
    // 	    hb.angle = 1.0;
    // 	    hb.energy = computeEnergy(distance, hb.angle);
	    
    // 	} else {
    // 	    std::vector<double> energies(4, 0.0);
    // 	    // std::cout << "Normal " << std::endl;
    // 	    k = buriedWaterVector.at(v1->info().id);
    // 	    l = buriedWaterVector.at(v2->info().id);
    // 	    OO = v1->point() - v2->point();
	    
    // 	    distance = 10.0*CGAL::sqrt(OO*OO);

    // 	    OH11 = *v1->info().h1 - v1->point();
    // 	    OH12 = *v1->info().h2 - v1->point();
    // 	    OH21 = *v2->info().h1 - v2->point();
    // 	    OH22 = *v2->info().h2 - v2->point();
	    
    // 	    OO = OO / CGAL::sqrt(OO*OO);
    // 	    OH11 = OH11 / CGAL::sqrt(OH11*OH11);
    // 	    OH12 = OH12 / CGAL::sqrt(OH12*OH12);
    // 	    OH21 = OH21 / CGAL::sqrt(OH21*OH21);
    // 	    OH22 = OH22 / CGAL::sqrt(OH22*OH22);

    // 	    std::vector<double> coss;
    // 	    coss.push_back(OH11 * OO);
    // 	    coss.push_back(OH12 * OO);
    // 	    coss.push_back(OH21 * -OO);
    // 	    coss.push_back(OH22 * -OO);
	    
    // 	    // std::cout << "cos11 : " << coss.at(0) << ", "
    // 	    // 	      << "cos12 : " << coss.at(1) << ", "
    // 	    // 	      << "cos21 : " << coss.at(2) << ", "
    // 	    // 	      << "cos22 : " << coss.at(3) << std::endl;
	    
    // 	    for (unsigned int i = 0; i < coss.size(); i++) {
    // 		if (coss.at(i) > 0.0) {
    // 		    energies.at(i) = computeEnergy(distance, coss.at(i));
    // 		    // std::cout << "energy" << i << " : " << energies.at(i) << std::endl;
    // 		}
    // 	    }
    // 	    auto result = std::max_element(energies.begin(), energies.end());

    // 	    int ind = std::distance(energies.begin(), result);
    // 	    if ( ind < 2 ) {
    // 		int temp;
    // 	        temp = k;
    // 		k = l;
    // 		l = temp;
    // 		way = true;
    // 	    }
    // 	    hb.length = distance;
    // 	    hb.angle = coss.at(ind);
    // 	    hb.energy = *result;
    // 	}

    // 	if ( hb.energy != 0.0 ) {
    // 	    // std::cout << "Distance : " << hb.length
    // 	    // 	      << ". Angle : " << hb.angle
    // 	    // 	      << ", Energy : " << hb.energy << std::endl;
    // 	    if (way) {
    // 		add_bidirectional_edge(v1->info().id, v2->info().id, hb, g, false);
     // 		add_bidirectional_edge(v1->info().id, v2->info().id, hb, gr, true);
    // 		plopVec.push_back(std::pair<std::pair<int, int>, double>(
    // 				      std::pair<int, int>(v1->info().id, v2->info().id),
    // 				      hb.energy));
    // 		hb_map.push_back(hb.energy);
    // 		hb.energy *= -1.0;
    // 		add_bidirectional_edge(v1->info().id, v2->info().id, hb, gp, true);
    // 	    } else {
    // 		add_bidirectional_edge(v2->info().id, v1->info().id, hb, g, false);
    // 		add_bidirectional_edge(v2->info().id, v1->info().id, hb, gr, true);
    // 		plopVec.push_back(std::pair<std::pair<int, int>, double>(
    // 				      std::pair<int, int>(v2->info().id, v1->info().id),
    // 				      hb.energy));
    // 		hb_map.push_back(hb.energy);
    // 		hb.energy *= -1.0;
    // 		add_bidirectional_edge(v1->info().id, v2->info().id, hb, gp, true);
    // 	    }
	    
    // 	    boost::add_edge(v1->info().id, v2->info().id, hb, gu);
    // 	    std::cout << hb.energy << ", ";
    // 	}
    // }

    // std::cout << std::endl;

    // std::cout << " >> Energies done " << std::endl;
    
    // /* Graph Analysis */
    // double flow = do_max_flow(g, Source_handle->info().id, Sink_handle->info().id);    
    // double flowf = do_max_flow(gr, Source_handle->info().id, Sink_handle->info().id);
    // double flowb = do_max_flow(gr, Sink_handle->info().id, Source_handle->info().id);
    
    // std::cout << " >> Flow done " << flow << " " << flowf << " " << flowb << std::endl;
    
    // /* Connected Component */
    // // std::vector<int> component (boost::num_vertices (g));
    // // size_t num_components = boost::connected_components (g, component.data()&component[0]);


    // // if (component.at(Source_handle->info().id) != component.at(Sink_handle->info().id)) {
    // // 	std::cout << " Number of Components : " << num_components 
    // // 		  << ", Component source : " << component.at(Source_handle->info().id)
    // // 		  << ", Component sink : " << component.at(Sink_handle->info().id) << std::endl;
    // // }
    

    // // int component_size = 0;
    // // int component_ind = component.at(Source_handle->info().id);
    // // for (int i = 0; i < component.size(); i++) {
    // // 	if (component.at(i) == component_ind) {
    // // 	    boost::add_vertex(Atom{i}, subg);
    // // 	    auto iter = boost::adjacent_vertices(i, g);
    // // 	    component_size++;
    // // 	}
    // // }
    
    // // int N = num_vertices(gr);
    // // std::vector<Graph::vertex_descriptor> pred(N, Graph::null_vertex()); 
    // // std::vector<int> dist(N, (std::numeric_limits < short >::max)());
    // // bool r = boost::bellman_ford_shortest_paths
    // // 	(gr, N, boost::weight_map(get(&HydrogenBond::energy, gr)).distance_map(dist.data()).
    // // 	 predecessor_map(pred.data()));


    // // for (auto &elem : dist) {
    // // 	std::cout << elem << " "; 
    // // }
    // // std::cout << "\n";
    
    // std::cout << " >> Components done " << component_size << std::endl;
    
    // int N = num_vertices(gu);
    // // Init pred vector 
    // std::vector<std::size_t> pred(N);
    // for (int i = 0; i< pred.size(); i++) {
    // 	pred.at(i) = i;
    // }

    // // Init distance vector with source distance to zero
    // std::vector<float> dist(N, (std::numeric_limits < short >::max)());
    // dist[Source_handle->info().id] = 0;
    // for (auto &elem : dist ) {
    // 	std::cout << elem << " ";
    // }    
    // std::cout << std::endl;
    
    // bool r = boost::bellman_ford_shortest_paths
    // 	(gp, N, boost::weight_map(get(&HydrogenBond::energy, gp)).distance_map(&dist[0]).
    // 	 predecessor_map(&pred[0]).root_vertex(Sink_handle->info().id));
    // if (r) {
    // 	for (auto &elem : pred ) {
    // 	    std::cout << elem << " ";
    // 	}
    // 	std::cout << std::endl;
    // 	for (auto &elem : dist ) {
    // 	    std::cout << elem << " ";
    // 	}    
    // 	std::cout << std::endl;
    // 	std::cout << " >> Path from " << Source_handle->info().id << " to "
    // 		  << Sink_handle->info().id << std::endl;
    // 	printPath(pred, Source_handle->info().id, Sink_handle->info().id);
    // 	std::cout << " >> Path done " << std::endl;
    // } else {
    // 	std::cout << " >> Path cycle problem " << std::endl;
    // }
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
    //}

filterData.startFrame(frnr, fr.time);
filterData.setPoint(0, results.numVertices);
filterData.setPoint(1, results.volume);
filterData.finishFrame();

graphData.startFrame(frnr, fr.time);
graphData.setPoint(0, results.numEdges);
graphData.setPoint(1, 1);
graphData.setPoint(2, 2);
graphData.finishFrame();

    // std::cout << " >> Output done" << std::endl; 
}


void WaterNetwork::finishAnalysis(int /*nframes*/)
{

}


void WaterNetwork::writeOutput()
{
    // outfilebinary.close();
    // outfilenormal.close();
}


int main(int argc, char *argv[])
{
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<WaterNetwork>(argc, argv);
}

