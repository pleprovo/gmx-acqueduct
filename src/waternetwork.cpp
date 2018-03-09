  

#include "waternetwork.hpp"
#include "gromacs/pbcutil/pbc.h"
#include <iostream>
#include <fstream>

#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>

#include <limits>


using Edge = std::pair<int, int>;

using HydrogenBond = std::pair<std::pair<int, int>, float>;

// Convert position from Gmx to CGAL type
template <class T>
std::vector<T> fromGmxtoCgalPosition(const gmx::ConstArrayRef<rvec> &coordinates,
				     const int increment=1) {

    std::vector<T> cgalPositionVector;   
    for (unsigned int i = 0; i < coordinates.size(); i += increment) {
	cgalPositionVector.push_back(T(coordinates.at(i)[XX],
				       coordinates.at(i)[YY],
				       coordinates.at(i)[ZZ]));
    }
    
    return cgalPositionVector;
}

struct edge_comparator {
    bool operator()(const std::pair<int, int> &a,
                    const std::pair<int, int> &b) const
	{
	    return less_comparator(std::minmax(a.first, a.second),
				   std::minmax(b.first, b.second));
	}

    std::less<std::pair<int, int> > less_comparator;
};

void checkPair(const gmx::AnalysisNeighborhoodPair &pair,
	       const gmx::Selection &waterSelection,
	       const gmx::ConstArrayRef<int> &mappedIds,
	       std::vector<Edge> &edgeVector,
	       const float &cutOffDistance = 0.35,
	       const float &cutOffAngle = 0.52) {
    rvec DH1, DH2;
    real angle1, angle2;
    real r_AD = norm(pair.dx());
    //if ( r_AD <= cutOffDistance ) {
	if (pair.refIndex() != pair.testIndex()) {
	    rvec_sub(waterSelection.position(mappedIds.at(pair.testIndex())+1).x(),
		     waterSelection.position(mappedIds.at(pair.testIndex())).x(),
		     DH1);
	    angle1 = gmx_angle(DH1, pair.dx());
	    if ( angle1 <= cutOffAngle ) {
		edgeVector.push_back(Edge(pair.refIndex(), pair.testIndex()));
	    }
	
	    rvec_sub(waterSelection.position(mappedIds.at(pair.testIndex())+2).x(),
		     waterSelection.position(mappedIds.at(pair.testIndex())).x(),
		     DH2);
	    angle2 = gmx_angle(DH2, pair.dx());
	    if ( angle2 <= cutOffAngle ) {
		edgeVector.push_back(Edge(pair.refIndex(), pair.testIndex()));
	    } 
	}
	//}
}

void checkHB(const gmx::AnalysisNeighborhoodPair &pair,
	       const gmx::Selection &waterSelection,
	       const gmx::ConstArrayRef<int> &mappedIds,
	       std::vector<HydrogenBond> &edgeVector,
	       const float &cutOffDistance = 0.35,
	       const float &cutOffAngle = 0.52) {
    rvec DH1, DH2;
    float energy_hb;
    /* Value for the hb potential that lead to average ~20kJ.mol-1*/
    float A = 0.0435; /* epsilon*sigma^6*sqrt(2/3) */
    float B = 0.5335; /* epsilon*sigma^4*sqrt(2/3) */
    real angle1, angle2;
    real r_AD = norm(pair.dx());
    HydrogenBond out;
    //if ( r_AD <= cutOffDistance ) {
	if (pair.refIndex() != pair.testIndex()) {
	    rvec_sub(waterSelection.position(mappedIds.at(pair.testIndex())+1).x(),
		     waterSelection.position(mappedIds.at(pair.testIndex())).x(),
		     DH1);
	    angle1 = gmx_angle(DH1, pair.dx());
	    if ( angle1 <= cutOffAngle ) {
		energy_hb = ((A/pow(r_AD, 6.0))-(B/pow(r_AD, 4.0)))*pow(cos(angle1), 4.0);
		out.first.first = pair.refIndex();
		out.first.second = pair.testIndex();
		out.second = energy_hb;
		edgeVector.push_back(out);
	    }
	
	    rvec_sub(waterSelection.position(mappedIds.at(pair.testIndex())+2).x(),
		     waterSelection.position(mappedIds.at(pair.testIndex())).x(),
		     DH2);
	    angle2 = gmx_angle(DH2, pair.dx());
	    if ( angle2 <= cutOffAngle ) {
		energy_hb = ((A/pow(r_AD, 6.0))-(B/pow(r_AD, 4.0)))*pow(cos(angle1), 4.0);
		out.first.first = pair.refIndex();
		out.first.second = pair.testIndex();
		out.second = energy_hb;
		edgeVector.push_back(out);
	    } 
	}
	//}
}

WaterNetwork::WaterNetwork()
    : TrajectoryAnalysisModule("waternetwork", "Water network analysis tool"),
      cutoff_(0.35)
{
    registerAnalysisDataset(&data_, "avedist");

    alphaShapeModule_ = std::make_shared<AlphaShapeModule>();
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

    options->addOption(gmx::SelectionOption("reference")
		       .store(&alphasel_).required()
		       .defaultSelectionText("Calpha")
		       .description("Reference group to calculate alpha shape (default C-alphas)"));
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
    watersel_.initOriginalIdsToGroup(top.topology(), INDEX_RES);

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
        plotm->setTitle("Average distance");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Distance (nm)");
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
    const gmx::Selection           &alphasel = pdata->parallelSelection(alphasel_);
    const gmx::Selection           &watersel = pdata->parallelSelection(watersel_);

    /* Get water Position and indices */
    gmx::ConstArrayRef<rvec> waterCoordinates = watersel.coordinates();

    /* Separate Oxygen coordinates and hydrogen coordinates */
    /* Get their indices */
    std::vector<rvec> oxygenVector(waterCoordinates.size()/3);
    std::vector<int> oxygenIndices, hydrogenIndices;
    
    for (unsigned int i = 0; i < oxygenVector.size(); i++) {
    	copy_rvec(waterCoordinates.at(3*i), oxygenVector.at(i));
    	oxygenIndices.push_back(3*i);
    	hydrogenIndices.push_back(3*i+1);
    	hydrogenIndices.push_back(3*i+2);
    }
    
    gmx::ConstArrayRef<rvec> oxygenCoordinates(oxygenVector);
    gmx::ConstArrayRef<int> oxygenIds(oxygenIndices);
    gmx::ConstArrayRef<int> hydrogenIds(hydrogenIndices);
      
    /* 
       In house Modules Initialization 
    */
    
    if (frnr == 0) {
    	this->dipoleModule_->initialise(waterCoordinates);
    	//this->lifetimeModule_->initialise(nb_frames, waterCoordinates.size()/3)
    }

    
    /* 
       Dipoles Stuff 
    */
    
    this->dipoleModule_->analyseFrame(waterCoordinates);

    
    /* 
       Alpha Shape stuff 
    */
    
    /* Convert the position into cgal Point */
    /* Compute the alpha shape for alpha bal of 1.0 nm */
    std::vector<int> buriedWaterVector;
    if (false) {
    	gmx::ConstArrayRef<rvec> alphaCoordinates = alphasel.coordinates();
    	std::vector<Point_3> alphaPoints = fromGmxtoCgalPosition<Point_3>(alphaCoordinates);
    	std::vector<Point_3> waterPoints = fromGmxtoCgalPosition<Point_3>(oxygenCoordinates);

    	alphaShapeModule_->build(alphaPoints, 1.0);    
    	buriedWaterVector = alphaShapeModule_->locate(waterPoints, TRUE);
	
    }

    
    /* 
       Graph Stuff 
    */
    
    /* Find all hydrogen bonds */
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
    while (pairSearch.findNextPair(&pair)) {
	checkPair(pair, watersel, oxygenIds, edgeVector);
	checkHB(pair, watersel, oxygenIds, HBVector);
    }

    std::set<int> sourceEdges;
    while (pairSearchSource.findNextPair(&pair)) {
	sourceEdges.insert(pair.refIndex());
    }
    
    std::set<int> sinkEdges;
    while (pairSearchSink.findNextPair(&pair)) {
	sinkEdges.insert(pair.refIndex());
    }
    
    float flow = 0.0;
    
    Graph g(oxygenIndices.size()+2);
    std::vector<edge_descriptor> reverseEdges;
    std::vector<float> capacity;
    
    for (auto hb : HBVector) {
	// edge_descriptor e1 = add_edge(hb.first.first, hb.first.second, g).first;
	// edge_descriptor e2 = add_edge(hb.first.second, hb.first.first, g).first;
	// put(boost::edge_capacity, g, e1, hb.second);
	// put(boost::edge_capacity, g, e2, hb.second);
	AddBidirectionalEdge(g, hb.first.first, hb.first.second, hb.second,
			     reverseEdges, capacity);
    }
    
    for (auto source : sourceEdges) {
	// edge_descriptor e1 = add_edge(oxygenIndices.size(), source, g).first;
	// edge_descriptor e2 = add_edge(source, s, g).first;
	// put(boost::edge_capacity, g, e1, std::numeric_limits<double>::infinity());
	// put(boost::edge_capacity, g, e2, std::numeric_limits<double>::infinity());
    }
	
    for (auto sink : sinkEdges) {
	// edge_descriptor e1 = add_edge(oxygenIndices.size()+1, sink, g).first;
	// edge_descriptor e2 = add_edge(sink, t, g).first;
	// put(boost::edge_capacity, g, e1, std::numeric_limits<double>::infinity());
	// put(boost::edge_capacity, g, e2, std::numeric_limits<double>::infinity());
    }
    int sourceId = 0;
    int sinkId = 1000;
    vertex_descriptor s = vertex(sourceId, g);
    vertex_descriptor t = vertex(sinkId, g);
    
    //flow = boykov_kolmogorov_max_flow(g, s, t);
    //Graph g(edgeVector.begin(), edgeVector.end(), oxygenCoordinates.size());
    //DiGraph dg(edgeVector.begin(), edgeVector.end(), oxygenCoordinates.size());
    
    /* Component analysis */
    //std::vector<ComponentGraph> components = connected_components_subgraphs(g);
    
    /* Motif and invariant search in connected components */

    /* Eigen value spectra analysis */
    
    /* Maximum Flow analysis between the sink and the source */
    
    /*
      Hydrogen Bonds Energy
    */
    float E = 0.0;
    for (auto hb: HBVector) {
	E += hb.second;
    }
    
    // Write graphml frame
    // Get all properties in dynamics properties
    
    // if ( frnr == 0) {
    // 	boost::dynamic_properties dp;
    // 	dp.property("vertex_index", get(boost::vertex_index_t(), dg));
    // 	dp.property("edge_index", get(boost::edge_index_t(), dg));
    // 	std::ofstream outFile("plop.dat");
    // 	std::stringstream oss;
    // 	write_graphml(oss, dg, dp, true);
    //     outFile << oss.str();
    // 	outFile.close();
    // }
    
    dh.startFrame(frnr, fr.time);
    dh.setPoint(0, 1.0*num_edges(g)/num_vertices(g));
    dh.setPoint(1, flow);
    dh.setPoint(2, 0.0);
    dh.finishFrame();
}


void
WaterNetwork::finishAnalysis(int /*nframes*/)
{

}


void
WaterNetwork::writeOutput()
{

}
