  

#include "waternetwork.hpp"
#include "gromacs/pbcutil/pbc.h"


using Edge = std::pair<int, int>;

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


WaterNetwork::WaterNetwork()
    : TrajectoryAnalysisModule("waternetwork", "Water network analysis tool"),
      cutoff_(0.25)
{
    registerAnalysisDataset(&data_, "avedist");

    alphaShapeModule_ = std::make_shared<AlphaShapeModule>();
    graphModule_      = std::make_shared<GraphModule>();
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

    options->addOption(gmx::DoubleOption("cutoff").store(&cutoff_)
		       .description("Cutoff for distance calculation (default 0.3)"));

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
    const gmx::Selection           &alphasel = pdata->parallelSelection(alphasel_);
    const gmx::Selection           &watersel = pdata->parallelSelection(watersel_);

    /* Get water Position and indices */
    gmx::ConstArrayRef<rvec> waterCoordinates = watersel.coordinates();
    gmx::ConstArrayRef<int> waterMappedIds = watersel.mappedIds();

    /* Separate Oxygen coordinates and hydrogen coordinates */
    /* Get their indices */
    std::vector<rvec> oxygenVector(waterCoordinates.size()/3);
    std::vector<rvec> hydrogenVector(2*(waterCoordinates.size()/3));
    std::vector<int> oxygenIndices, hydrogenIndices;
    
    for (unsigned int i = 0; i < oxygenVector.size(); i++) {
    	copy_rvec(waterCoordinates.at(3*i), oxygenVector.at(i));
    	copy_rvec(waterCoordinates.at(3*i+1), hydrogenVector.at(2*i));
    	copy_rvec(waterCoordinates.at(3*i+2), hydrogenVector.at(2*i+1));
    	oxygenIndices.push_back(3*i);
    	hydrogenIndices.push_back(3*i+1);
    	hydrogenIndices.push_back(3*i+2);
    }
    
    gmx::ConstArrayRef<rvec> oxygenCoordinates(oxygenVector);
    gmx::ConstArrayRef<rvec> hydrogenCoordinates(hydrogenVector);
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
    if (true) {
    	gmx::ConstArrayRef<rvec> alphaCoordinates = alphasel.coordinates();
    	std::vector<Point_3> alphaPoints = fromGmxtoCgalPosition<Point_3>(alphaCoordinates);
    	std::vector<Point_3> waterPoints = fromGmxtoCgalPosition<Point_3>(oxygenCoordinates);

    	alphaShapeModule_->build(alphaPoints, 1.0);    
    	buriedWaterVector = alphaShapeModule_->locate(waterPoints, TRUE);
	
    }

    
    /* 
       Graph Stuff 
    */
    
    /* Build set of H--O pairs */
    std::set<Edge, edge_comparator> edgeSet;
    gmx::AnalysisNeighborhoodPositions waterPos(watersel);
    /*
    gmx::AnalysisNeighborhoodPositions oxygenPos(oxygenCoordinates.data(),
    						 oxygenCoordinates.size());
    gmx::AnalysisNeighborhoodPositions hydrogenPos(hydrogenCoordinates.data(),
    						   hydrogenCoordinates.size());
    */
    gmx::AnalysisNeighborhoodPositions oxygenPos = waterPos.indexed(oxygenIds);
    gmx::AnalysisNeighborhoodPositions hydrogenPos = waterPos.indexed(hydrogenIds);
    gmx::AnalysisNeighborhoodSearch search = nb_.initSearch(pbc, hydrogenPos);
    gmx::AnalysisNeighborhoodPairSearch pairSearch = search.startPairSearch(oxygenPos);
    gmx::AnalysisNeighborhoodPair pair;

    while (pairSearch.findNextPair(&pair))
    {
    	if (norm(pair.dx()) > 0.12) {
    	    edgeSet.insert(Edge(pair.refIndex(), pair.testIndex()));
    	}
    }
    
    std::vector<Edge> edgeVector(edgeSet.begin(), edgeSet.end());

    Graph g(edgeVector.begin(), edgeVector.end(), oxygenCoordinates.size());
    
    /* Get connected components and compute their dipoles moment orientation */
    
    /* Size distribution in bulk and buried waters */
    
    /* Dipoles distribution in bulk and buried waters */

    /* Component analysis */
    std::vector<ComponentGraph> components = connected_components_subgraphs(g);
    
    /* Motif and invariant search in connected components */


    /* Eigen value spectra analysis */
    
    
    /* 
       Data Point Storage Stuff 
    */
    
    dh.startFrame(frnr, fr.time);
    dh.setPoint(0, buriedWaterVector.size());
    dh.setPoint(1, edgeSet.size());
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
