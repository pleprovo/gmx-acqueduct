  
// #include "cgal.hpp"
#include "alphashape.hpp"

#include <iostream>
#include <fstream>
#include <algorithm> 

// Convert position from Gmx to CGAL type
template <class T>
std::vector<T> fromGmxtoCgalPosition(const gmx::ArrayRef<const rvec> &coordinates,
				     const int increment=1)
{
    std::vector<T> cgalPositionVector;   
    for (unsigned int i = 0; i < coordinates.size(); i += increment)
    {
	cgalPositionVector.push_back(T(coordinates.at(i)[XX],
				       coordinates.at(i)[YY],
				       coordinates.at(i)[ZZ]));
    }
    return cgalPositionVector;
}


AlphaShape::AlphaShape() : lifetimeModule_(new gmx::AnalysisDataLifetimeModule())
{
    alphaValue_ = 1.0;
	
    numFrameValue_ = 250;
    registerAnalysisDataset(&waterData_, "aveWater");
    registerAnalysisDataset(&volumeData_, "aveVolume");
    
    lifetimeData_.addModule(lifetimeModule_);
    
    registerAnalysisDataset(&lifetimeData_, "data");
    registerBasicDataset(lifetimeModule_.get(), "lifetime"); 
}


void AlphaShape::initOptions(gmx::IOptionsContainer          *options,
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

    settings->setHelpText(desc);

    options->addOption(gmx::FileNameOption("o").filetype(gmx::eftPlot).outputFile()
		       .store(&filenameWater_).defaultBasename("buried")
		       .description("Number of buried water as function of time"));
    options->addOption(gmx::FileNameOption("ov").filetype(gmx::eftPlot).outputFile()
	               .store(&filenameVolume_).defaultBasename("volume")
		       .description("Write the OFF file of the last alpha shape"));
    options->addOption(gmx::FileNameOption("ol").filetype(gmx::eftPlot).outputFile()
	               .store(&filenameLifetime_).defaultBasename("life")
		       .description("Lifetime of water molecules within the surface"));  
    options->addOption(gmx::SelectionOption("Alpha Points")
		       .storeVector(&selectionListAlpha_)
		       .required().multiValue()
		       .description("Reference group to calculate alpha shape)"));
    options->addOption(gmx::SelectionOption("Waters").store(&selectionWater_)
		       .required().defaultSelectionText("Water")
		       .description("Groups to calculate graph properties (default Water)"));
    
    options->addOption(gmx::DoubleOption("alpha").store(&alphaValue_)
		       .description("Alpha Value for the Alpha Shape computation"));
    options->addOption(gmx::DoubleOption("nf").store(&numFrameValue_)
		       .description("Number fof frames to consider for lifetime"));
    
}

void AlphaShape::optionsFinished(gmx::TrajectoryAnalysisSettings *settings)
{

    settings->setFlag(gmx::TrajectoryAnalysisSettings::efRequireTop);
    settings->setFlag(gmx::TrajectoryAnalysisSettings::efUseTopX);

}

void AlphaShape::initAnalysis(const gmx::TrajectoryAnalysisSettings &settings,
			      const gmx::TopologyInformation        &top)
{
    as_ = std::make_shared<AlphaShapeSurface>();
    as_->setAlphaValue(alphaValue_);

/* Set the number of column to store time dependent data */
    waterData_.setColumnCount(0, selectionListAlpha_.size());
    volumeData_.setColumnCount(0, selectionListAlpha_.size());

    // alphaShapeModulePtr_->setAlpha(alphaValue_);
    std::cout << "Using Alpha Value : " << alphaValue_ << std::endl;
    std::cout << "Using persistance value : " << numFrameValue_ << std::endl;

    
    /* Init the Plot module for the time dependent data */
    if (!filenameWater_.empty())
    {
	gmx::AnalysisDataPlotModulePointer plotm(
	    new gmx::AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(filenameWater_);
        plotm->setTitle("Buried Water");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Number of Water");
	for (size_t i = 0; i < selectionListAlpha_.size(); i++)
	{
	    plotm->appendLegend(std::string(selectionListAlpha_.at(i).name()));
	}
        waterData_.addModule(plotm);
    }

    if (!filenameVolume_.empty())
    {
	gmx::AnalysisDataPlotModulePointer plotm(
	    new gmx::AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(filenameVolume_);
        plotm->setTitle("Volume");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Volume (nm^3)");
	for (size_t i = 0; i < selectionListAlpha_.size(); i++)
	{
	    plotm->appendLegend(std::string(selectionListAlpha_.at(i).name()));
	}
        volumeData_.addModule(plotm);
    }
    
    lifetimeData_.setDataSetCount(selectionListAlpha_.size());
    for (size_t g = 0; g < selectionListAlpha_.size(); g++)
    {
    	lifetimeData_.setColumnCount(g, selectionWater_.posCount()/3);
    }
    lifetimeModule_->setCumulative(true);
    if (!filenameLifetime_.empty())
    {
    	gmx::AnalysisDataPlotModulePointer plot(
    	    new gmx::AnalysisDataPlotModule(settings.plotSettings()));
    	plot->setFileName(filenameLifetime_);
    	plot->setTitle("Lifetime of water in the surface");
    	plot->setXAxisIsTime();
    	plot->setYLabel("Occupancy");

    	for (size_t g = 0; g < selectionListAlpha_.size(); g++)
    	{
    	    plot->appendLegend(std::string(selectionListAlpha_.at(g).name()));
    	}
    	lifetimeModule_->addModule(plot);
	
    }
}


void AlphaShape::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
			   gmx::TrajectoryAnalysisModuleData *pdata)
{
    gmx::AnalysisDataHandle waterDataHandle = pdata->dataHandle(waterData_);
    gmx::AnalysisDataHandle volumeDataHandle = pdata->dataHandle(volumeData_);
    gmx::AnalysisDataHandle lifetimeDataHandle = pdata->dataHandle(lifetimeData_);
    
    const gmx::SelectionList &selectionListAlpha = pdata->parallelSelections(selectionListAlpha_);
    const gmx::Selection &selectionWater = pdata->parallelSelection(selectionWater_);

    /* Get water Position and indices */
    const gmx::ArrayRef<const rvec>& waterCoordinates = selectionWater.coordinates();

    /* Convert gromacs rvec position to cgal Point class */
    std::vector<Point> oxygens = fromGmxtoCgalPosition<Point>(waterCoordinates, 3);
    
    waterDataHandle.startFrame(frnr, fr.time);
    volumeDataHandle.startFrame(frnr, fr.time);
    lifetimeDataHandle.startFrame(frnr, fr.time);
    for (unsigned int i = 0; i < selectionListAlpha.size(); i++)
    {	
    	std::vector<Point> alphaPoints = fromGmxtoCgalPosition<Point>(selectionListAlpha.at(i).coordinates());	

    	std::vector<int> selected;
	as_->make(alphaPoints);
	int num_filtered = as_->locate(oxygens, selected);

    	waterDataHandle.setPoint(i, num_filtered);
	
	volumeDataHandle.setPoint(i, as_->volume());

       	lifetimeDataHandle.selectDataSet(i);
    	for (unsigned int j = 0; j < oxygens.size(); ++j)
    	{
    	    lifetimeDataHandle.setPoint(j, 0);
    	}
    	for (int id : selected)
    	{
    	    lifetimeDataHandle.setPoint(id, 1);
    	}
	
    }
    lifetimeDataHandle.finishFrame();
    waterDataHandle.finishFrame();
    volumeDataHandle.finishFrame();
}


void AlphaShape::finishAnalysis(int nframes)
{

}


void AlphaShape::writeOutput()
{

}

int main(int argc, char *argv[])
{
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<AlphaShape>(argc, argv);
}

