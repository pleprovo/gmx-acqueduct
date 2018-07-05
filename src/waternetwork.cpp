  

#include "waternetwork.hpp"
#include "gromacs/pbcutil/pbc.h"

#include <iostream>
#include <fstream>
#include <limits>


using Edge = std::pair<unsigned int, unsigned int>;

using HydrogenBond = std::pair<Edge, float>;

float energy_hb(const rvec& donor, const rvec& acceptor, const rvec& hydrogen) {
    return 0.0;
}

void checkHB(const gmx::AnalysisNeighborhoodPair &pair,
	     const gmx::Selection &waterSelection,
	     const gmx::ConstArrayRef<int> &mappedIds,
	     std::vector<HydrogenBond> &edgeVector,
 	     const float &cutOffAngle = 0.52)
{
    rvec DH1, DH2;
    float energy_hb;
    float C = 3855; /* epsilon*sigma^6*sqrt(2/3) */
    float D = 738; /* epsilon*sigma^4*sqrt(2/3) */
    real angle1, angle2;
    real r_AD = 10.0*norm(pair.dx());
    HydrogenBond out;
    if (pair.refIndex() != pair.testIndex()) {
	rvec_sub(waterSelection.position(mappedIds.at(pair.testIndex())+1).x(),
		 waterSelection.position(mappedIds.at(pair.testIndex())).x(),
		 DH1);
	angle1 = gmx_angle(DH1, pair.dx());
	
	//if ( angle1 > 0.0 ) {
	energy_hb = -((C/pow(r_AD, 6.0))-(D/pow(r_AD, 4.0))*pow(cos(angle1), 4.0));
	out.first.first = pair.refIndex();
	out.first.second = pair.testIndex();
	out.second = energy_hb;
	edgeVector.push_back(out);
	    //}
	
	rvec_sub(waterSelection.position(mappedIds.at(pair.testIndex())+2).x(),
		 waterSelection.position(mappedIds.at(pair.testIndex())).x(),
		 DH2);
	angle2 = gmx_angle(DH2, pair.dx());
	
	//if ( angle2 > 0.0 ) {
	energy_hb = -((C/pow(r_AD, 6.0))-(D/pow(r_AD, 4.0))*pow(cos(angle2), 4.0));
	out.first.first = pair.refIndex();
	out.first.second = pair.testIndex();
	out.second = energy_hb;
	edgeVector.push_back(out);
	    //} 
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
    
    //int nextEdgeId = num_edges(graph);
    int nextEdgeId = capacity.size();
    edge_descriptor edge;
    bool inserted;

    boost::tie(edge,inserted) = add_edge(source, target, nextEdgeId, graph);
    if(!inserted) {
        std::cerr << "Not inserted!" << std::endl;
    }
    edge_descriptor reverseEdge = add_edge(target, source, nextEdgeId + 1, graph).first;
    reverseEdges.push_back(reverseEdge);
    reverseEdges.push_back(edge);
    capacity.push_back(weight);

    // Not sure what to do about reverse edge weights
    capacity.push_back(weight);
    // capacity.push_back(0);
}

WaterNetwork::WaterNetwork()
    : TrajectoryAnalysisModule("waternetwork", "Water network analysis tool"),
      cutoff_(0.5)
{
    registerAnalysisDataset(&data_, "avedist");

    graphModulePtr_ = std::make_shared<GraphModule>();
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

    options->addOption(gmx::FileNameOption("o")
		       .filetype(gmx::eftPlot).outputFile()
		       .store(&fnDist_).defaultBasename("avedist")
		       .description("Collection of analysis properties through time"));

    options->addOption(gmx::SelectionOption("select")
		       .store(&solvent_).required()
		       .defaultSelectionText("Water")
		       .description("Groups to calculate graph properties (default Water)"));
    
    options->addOption(gmx::SelectionOption("alpha")
		       .store(&calpha_).required()
		       .defaultSelectionText("Protein")
		       .description(""));
    
    // Source and sink option
    options->addOption(gmx::SelectionOption("source")
		       .store(&source_)
		       .description("Define a group as a source for maximum flow analysis, must be used with -sink option"));

    options->addOption(gmx::SelectionOption("sink")
		       .store(&sink_)
		       .description("Define a group as a sink for maximum flow analysis, must be used with -source option"));
    
    options->addOption(gmx::DoubleOption("cutoff").store(&cutoff_)
		       .description("Cutoff for distance calculation (default 0.25 nm)"));

    settings->setFlag(gmx::TrajectoryAnalysisSettings::efRequireTop);
    settings->setFlag(gmx::TrajectoryAnalysisSettings::efUseTopX);
}


void WaterNetwork::initAnalysis(const gmx::TrajectoryAnalysisSettings &settings,
			   const gmx::TopologyInformation        &top)
{
    /* Init Selection */

    solvent_.initOriginalIdsToGroup(top.topology(), INDEX_RES);
    nb_water_ = solvent_.posCount()/3;
    
    /* Init neihborsearch cutoff value */
    nb1_.setCutoff(cutoff_);
    nb2_.setCutoff(0.30);
    nb1_.setMode(gmx::AnalysisNeighborhood::SearchMode::eSearchMode_Grid);
    nb2_.setMode(gmx::AnalysisNeighborhood::SearchMode::eSearchMode_Grid);
    
    /* Create ConstArray of indices for oxygen atoms */
    oxygenIndices_.reserve(nb_water_);
    for (int i = 0; i < nb_water_; i++) {
    	oxygenIndices_.push_back(3*i);
    }

    std::cout << "THE LAST WATER : " << oxygenIndices_.back() << std::endl;
    std::cout << "NUMBER OF WATER : " << oxygenIndices_.size() << std::endl;
    
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

    gmx::ConstArrayRef<int> oxygenArrayRef(oxygenIndices_);
    
    /* Find all hydrogen bonds */
    gmx::AnalysisNeighborhoodPair pair;
    gmx::AnalysisNeighborhoodPositions sourcePos(sourcesel);
    gmx::AnalysisNeighborhoodPositions calphaPos(calphasel);
    gmx::AnalysisNeighborhoodPositions sinkPos(sinksel);    
    gmx::AnalysisNeighborhoodPositions waterPos(watersel);
    
    gmx::AnalysisNeighborhoodPositions oxygenPos = waterPos.indexed(oxygenArrayRef);
    
    gmx::AnalysisNeighborhoodSearch searchAlpha = nb1_.initSearch(pbc, oxygenPos);
    gmx::AnalysisNeighborhoodPairSearch pairSearchAlpha = searchAlpha.startPairSearch(calphaPos);
    std::set<int> neighborSet;
    std::vector<int> neighborVec;
    while (pairSearchAlpha.findNextPair(&pair)) {
        neighborSet.insert(oxygenArrayRef.at(pair.refIndex()));
    }
    // std::cout << waterSet.size() << std::endl; 
    std::copy(neighborSet.begin(), neighborSet.end(), std::back_inserter(neighborVec));
    gmx::ConstArrayRef<int> neighborArrayRef(neighborVec);
    
    gmx::AnalysisNeighborhoodPositions neighborPos = waterPos.indexed(neighborArrayRef);
    gmx::AnalysisNeighborhoodSearch search1 = nb1_.initSearch(pbc, neighborPos);
    gmx::AnalysisNeighborhoodPairSearch pairSearch = search1.startPairSearch(neighborPos);
	
    gmx::AnalysisNeighborhoodSearch search2 = nb2_.initSearch(pbc, neighborPos);    
    gmx::AnalysisNeighborhoodPairSearch pairSearchSource = search2.startPairSearch(sourcePos);
    gmx::AnalysisNeighborhoodPairSearch pairSearchSink = search2.startPairSearch(sinkPos);
    
    

    std::set<int> sourceEdges;
    std::set<int> sinkEdges;
    std::vector<HydrogenBond> HBVector;
    
    while (pairSearch.findNextPair(&pair)) {
    	checkHB(pair, watersel, neighborArrayRef, HBVector);
    }

    while (pairSearchSource.findNextPair(&pair)) {
    	sourceEdges.insert(pair.refIndex());
    }

    while (pairSearchSink.findNextPair(&pair)) {
    	sinkEdges.insert(pair.refIndex());
    }

    float flow = 0.0;
    int sourceId = nb_water_;
    int sinkId = nb_water_+1;
    Graph g(nb_water_+2);
    
    std::vector<edge_descriptor> reverseEdges;
    std::vector<float> capacity;
    for (const auto& hb : HBVector) {
    	AddBidirectionalEdge(g, hb.first.first, hb.first.second, hb.second,
    			     reverseEdges, capacity);	
    }
    
    EdgeProperties infinityEdge;
    infinityEdge.capacity = 1000.0;
    for (const auto& source : sourceEdges) {
    	AddBidirectionalEdge(g, sourceId, source, 1000.0, reverseEdges, capacity);
    }
    
    for (const auto& sink : sinkEdges) {
    	AddBidirectionalEdge(g, sinkId, sink, 1000.0, reverseEdges, capacity);
    }
    
    vertex_descriptor s = vertex(sourceId, g);
    vertex_descriptor t = vertex(sinkId, g);    
    std::vector<float> residual_capacity(num_edges(g), 0.0);
    
    flow = boost::boykov_kolmogorov_max_flow(g,
    	   boost::make_iterator_property_map(&capacity[0], get(boost::edge_index, g)),
    	   boost::make_iterator_property_map(&residual_capacity[0], get(boost::edge_index, g)),
    	   boost::make_iterator_property_map(&reverseEdges[0], get(boost::edge_index, g)),
    					     get(boost::vertex_index, g), s, t);

    
    dh.startFrame(frnr, fr.time);
    dh.setPoint(0, sourceEdges.size());
    dh.setPoint(1, sinkEdges.size());
    dh.setPoint(2, 0 /*waterVectorRef.size()*/);
    dh.setPoint(3, flow);
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

