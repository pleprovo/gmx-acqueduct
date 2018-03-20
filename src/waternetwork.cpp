  

#include "waternetwork.hpp"
#include "gromacs/pbcutil/pbc.h"
#include <iostream>
#include <fstream>

#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>

#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <limits>

#include <ctime>

using Edge = std::pair<unsigned int, unsigned int>;

using HydrogenBond = std::pair<Edge, float>;

void checkHB(const gmx::AnalysisNeighborhoodPair &pair,
	     const gmx::Selection &waterSelection,
	     const gmx::ConstArrayRef<int> &mappedIds,
	     std::vector<HydrogenBond> &edgeVector,
	     const float &cutOffAngle = 0.52)
{
    rvec DH1, DH2;
    float energy_hb;
    /* Value for the hb potential that lead to average ~20kJ.mol-1*/
    float C = 3855; /* epsilon*sigma^6*sqrt(2/3) */
    float D = 738; /* epsilon*sigma^4*sqrt(2/3) */
    real angle1, angle2;
    real r_AD = 10.0*norm(pair.dx());
    HydrogenBond out;
    if (pair.refIndex() != pair.testIndex())
    {
	rvec_sub(waterSelection.position(mappedIds.at(pair.testIndex())+1).x(),
		 waterSelection.position(mappedIds.at(pair.testIndex())).x(),
		 DH1);
	angle1 = gmx_angle(DH1, pair.dx());
	
	if ( angle1 <= cutOffAngle )
	{
	    energy_hb = -((C/pow(r_AD, 6.0))-(D/pow(r_AD, 4.0)))*pow(cos(angle1), 4.0);
	    out.first.first = pair.refIndex();
	    out.first.second = pair.testIndex();
	    out.second = energy_hb;
	    edgeVector.push_back(out);
	}
	
	rvec_sub(waterSelection.position(mappedIds.at(pair.testIndex())+2).x(),
		 waterSelection.position(mappedIds.at(pair.testIndex())).x(),
		 DH2);
	angle2 = gmx_angle(DH2, pair.dx());
	
	if ( angle2 <= cutOffAngle )
	{
	    energy_hb = -((C/pow(r_AD, 6.0))-(D/pow(r_AD, 4.0)))*pow(cos(angle2), 4.0);
	    out.first.first = pair.refIndex();
	    out.first.second = pair.testIndex();
	    out.second = energy_hb;
	    edgeVector.push_back(out);
	} 
    }
}

void AddBidirectionalEdge(Graph &graph,
			  const unsigned int &source,
			  const unsigned int &target,
			  const float &weight,
                          std::vector<edge_descriptor>& reverseEdges,
			  std::vector<float>& capacity)
{
    // Add edges between grid vertices. We have to create the edge and the reverse edge,
    // then add the reverseEdge as the corresponding reverse edge to 'edge', and then add 'edge'
    // as the corresponding reverse edge to 'reverseEdge'
    int nextEdgeId = num_edges(graph);

    edge_descriptor edge;
    bool inserted;

    boost::tie(edge,inserted) = add_edge(source, target, nextEdgeId, graph);
    if(!inserted)
    {
        std::cerr << "Not inserted!" << std::endl;
    }
    edge_descriptor reverseEdge = add_edge(target, source, nextEdgeId + 1, graph).first;
    reverseEdges.push_back(reverseEdge);
    reverseEdges.push_back(edge);
    capacity.push_back(weight);

    // Not sure what to do about reverse edge weights
    capacity.push_back(weight);
//    capacity.push_back(0);
}


WaterNetwork::WaterNetwork()
    : TrajectoryAnalysisModule("waternetwork", "Water network analysis tool"),
      cutoff_(0.35)
{
    registerAnalysisDataset(&data_, "avedist");

    dipoleModule_     = std::make_shared<DipoleModule>();
    lifetimeModule_   = std::make_shared<LifetimeModule>();
}

void
WaterNetwork::initOptions(gmx::Options                    *options,
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

    options->addOption(gmx::FileNameOption("o")
		       .filetype(gmx::eftPlot).outputFile()
		       .store(&fnDist_).defaultBasename("avedist")
		       .description("Collection of analysis properties through time"));

    options->addOption(gmx::SelectionOption("select")
		       .store(&watersel_).required()
		       .defaultSelectionText("Water")
		       .description("Groups to calculate graph properties (default Water)"));
    
    // Source and sink option
    options->addOption(gmx::SelectionOption("source")
		       .store(&sourcesel_)
		       .description("Define a group as a source for maximum flow analysis, must be used with -sink option"));

    options->addOption(gmx::SelectionOption("sink")
		       .store(&sinksel_)
		       .description("Define a group as a sink for maximum flow analysis, must be used with -source option"));
    
    options->addOption(gmx::DoubleOption("cutoff").store(&cutoff_)
		       .description("Cutoff for distance calculation (default 0.25 nm)"));

    settings->setFlag(gmx::TrajectoryAnalysisSettings::efRequireTop);
    settings->setFlag(gmx::TrajectoryAnalysisSettings::efUseTopX);
}


void
WaterNetwork::initAnalysis(const gmx::TrajectoryAnalysisSettings &settings,
			   const gmx::TopologyInformation        &top)
{
    /* Init Selection */
    int nb_water = watersel_.initOriginalIdsToGroup(top.topology(), INDEX_RES);
    std::cout << nb_water << std::endl;
    /* Init neihborsearch cutoff value */
    nb_.setCutoff(cutoff_);

    /* Set the number of column to store time dependent data */
    data_.setColumnCount(0, 3);

    /* Init the average module  */ 
    avem_.reset(new gmx::AnalysisDataAverageModule());
    data_.addModule(avem_);

    /* Init the Plot module for the time dependent data */
    if (!fnDist_.empty())
    {
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
    const gmx::Selection           &sourcesel = pdata->parallelSelection(sourcesel_);
    const gmx::Selection           &sinksel = pdata->parallelSelection(sinksel_);
    const gmx::Selection           &watersel = pdata->parallelSelection(watersel_);

    /* Separate Oxygen coordinates and hydrogen coordinates */
    /* Get their indices */
    std::vector<int> oxygenIndices, hydrogenIndices;
    
    for (unsigned int i = 0; i < watersel.posCount()/3; i++)
    {
    	oxygenIndices.push_back(3*i);
    	hydrogenIndices.push_back(3*i+1);
    	hydrogenIndices.push_back(3*i+2);
    }
    
    gmx::ConstArrayRef<int> oxygenIds(oxygenIndices);
    gmx::ConstArrayRef<int> hydrogenIds(hydrogenIndices);
      
    /* 
       In house Modules Initialization 
    */
    
    if (frnr == 0)
    {
    	// this->dipoleModule_->initialise(waterCoordinates);
    	// this->lifetimeModule_->initialise(nb_frames, waterCoordinates.size()/3)
    }

    
    /* 
       Dipoles Stuff 
    */
    
    // this->dipoleModule_->analyseFrame(waterCoordinates);

    
    /* 
       Graph Stuff 
    */
    
    /* Find all hydrogen bonds */

    clock_t begin = clock();
    
    gmx::AnalysisNeighborhoodPositions sourcePos(sourcesel);
    gmx::AnalysisNeighborhoodPositions sinkPos(sinksel);
    gmx::AnalysisNeighborhoodPositions waterPos(watersel);
    gmx::AnalysisNeighborhoodPositions oxygenPos = waterPos.indexed(oxygenIds);
    
    gmx::AnalysisNeighborhoodSearch search = nb_.initSearch(pbc, oxygenPos);    
    gmx::AnalysisNeighborhoodPairSearch pairSearchSource = search.startPairSearch(sourcePos);
    gmx::AnalysisNeighborhoodPairSearch pairSearchSink = search.startPairSearch(sinkPos);
    gmx::AnalysisNeighborhoodPairSearch pairSearch = search.startPairSearch(oxygenPos);
    
    gmx::AnalysisNeighborhoodPair pair;
    
    std::vector<Edge> edgeVector;
    std::vector<HydrogenBond> HBVector;
    while (pairSearch.findNextPair(&pair))
    {
	checkHB(pair, watersel, oxygenIds, HBVector);
	//megaHB(pair, watersel, oxygenIds, g, reverseEdges, capacity);
    }

    std::set<int> sourceEdges;
    while (pairSearchSource.findNextPair(&pair))
    {
	sourceEdges.insert(pair.refIndex());
    }
    
    std::set<int> sinkEdges;
    while (pairSearchSink.findNextPair(&pair))
    {
	sinkEdges.insert(pair.refIndex());
    }
      
    clock_t top_1 = clock();
    
    float flow = 0.0;
    int sourceId = 0;
    int sinkId = 1;
    Graph g(oxygenIndices.size()+2);
    std::vector<edge_descriptor> reverseEdges;
    std::vector<float> capacity;
    
    for (const auto& hb : HBVector)
    {
    	AddBidirectionalEdge(g, hb.first.first, hb.first.second, hb.second,
    			     reverseEdges, capacity);
    }
    
    for (const auto& source : sourceEdges)
    {
    	AddBidirectionalEdge(g, sourceId, source, 1000.0 /*std::numeric_limits<float>::infinity()*/,
    			     reverseEdges, capacity);
    }
	
    for (const auto& sink : sinkEdges)
    {
    	AddBidirectionalEdge(g, sinkId, sink, 1000.0 /*std::numeric_limits<float>::infinity()*/,
    			     reverseEdges, capacity);
    }

    // vertex_descriptor s = vertex(sourceId, g);
    // vertex_descriptor t = vertex(sinkId, g);
    // std::vector<float> residual_capacity(num_edges(g), 0);

    // flow = boost::boykov_kolmogorov_max_flow(g,
    // 	   boost::make_iterator_property_map(&capacity[0], get(boost::edge_index, g)),
    // 	   boost::make_iterator_property_map(&residual_capacity[0], get(boost::edge_index, g)),
    // 	   boost::make_iterator_property_map(&reverseEdges[0], get(boost::edge_index, g)),
    // 	   get(boost::vertex_index, g), s, t);
  
    // // std::vector<vertex_descriptor> p(num_vertices(g));
    // // boost::dijkstra_shortest_paths(g, s, boost::predecessor_map(&p[0]));
    
    clock_t top_2 = clock();
    
    float E = 0.0;
    float maxE = 0.0;
    float minE = 0.0;
    for (auto hb: HBVector)
    {
    	if ( hb.second > maxE )
    	{
    	    maxE = hb.second;
    	}
    	if ( hb.second < minE )
    	{
    	    minE = hb.second;
    	}
    	E += hb.second;
    }

    
    
    this->search_time_ += double(top_1 - begin) / CLOCKS_PER_SEC;
    this->graph_time_ += double(top_2 - top_1) / CLOCKS_PER_SEC;
    
    dh.startFrame(frnr, fr.time);
    dh.setPoint(0, sourceEdges.size());
    dh.setPoint(1, sinkEdges.size());
    dh.setPoint(2, flow);
    dh.finishFrame();
    
}


void
WaterNetwork::finishAnalysis(int /*nframes*/)
{
    double total = this->search_time_ + this->graph_time_ + this-> analysis_time_;
    std::cout << 100*this->search_time_/total << " / "
	      << 100*this->graph_time_/total << " / "
	      << 100*this->analysis_time_/total << std::endl;
}


void
WaterNetwork::writeOutput()
{

}

int
main(int argc, char *argv[])
{
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<WaterNetwork>(argc, argv);
}
