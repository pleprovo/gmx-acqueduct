  

#include "alphashape.hpp"

#include <iostream>
#include <fstream>
#include <algorithm> 

// Convert position from Gmx to CGAL type
template <class T>
std::vector<T> fromGmxtoCgalPosition(const gmx::ConstArrayRef<rvec> &coordinates,
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

template <typename T>
std::vector<std::vector<T> > transpose_2D(const std::vector<std::vector<T> > &mat)
{
    std::vector<std::vector<T> > out(mat.at(0).size(), std::vector<T>(mat.size()));
    for (unsigned int i = 0; i < mat.size(); ++i)
    {
	for (unsigned int j = 0; j < mat.at(0).size(); ++j)
	{
	    out.at(j).at(i) = mat.at(i).at(j);
	}
    }
    return out;
}

AlphaShape::AlphaShape()
    : TrajectoryAnalysisModule("AlphaShape", "Protein Surface analysis")
{
    alphaValue_ = 1.0;
    numFrameValue_ = 250;
    registerAnalysisDataset(&waterData_, "aveWater");
    registerAnalysisDataset(&volumeData_, "aveVolume");
    alphaShapeModulePtr_ = std::make_shared<AlphaShapeModule>();
}


void AlphaShape::initOptions(gmx::Options                    *options,
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
		       .store(&filenameWater_).defaultBasename("buried")
		       .description("Number of buried water as function of time"));
    options->addOption(gmx::FileNameOption("ov").filetype(gmx::eftPlot).outputFile()
	               .store(&filenameVolume_).defaultBasename("volume")
		       .description("Write the OFF file of the last alpha shape"));
    options->addOption(gmx::FileNameOption("os").filetype(gmx::eftUnknown).outputFile()
	               .store(&filenameStats_).defaultBasename("stats")
		       .description("Write the stats of buried water and most frequent water"));  
    options->addOption(gmx::SelectionOption("Alpha Points").storeVector(&selectionListAlpha_)
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


void AlphaShape::initAnalysis(const gmx::TrajectoryAnalysisSettings &settings,
			 const gmx::TopologyInformation        &top)
{
    /* Set the number of column to store time dependent data */
    waterData_.setColumnCount(0, selectionListAlpha_.size());
    volumeData_.setColumnCount(0, selectionListAlpha_.size());

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
    
}


void AlphaShape::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
			   gmx::TrajectoryAnalysisModuleData *pdata)
{
    gmx::AnalysisDataHandle waterDataHandle = pdata->dataHandle(waterData_);
    gmx::AnalysisDataHandle volumeDataHandle = pdata->dataHandle(volumeData_);
    
    const gmx::SelectionList &selectionListAlpha = pdata->parallelSelections(selectionListAlpha_);
    const gmx::Selection &selectionWater = pdata->parallelSelection(selectionWater_);

    /* Get water Position and indices */
    gmx::ConstArrayRef<rvec> waterCoordinates = selectionWater.coordinates();

    /* Convert gromacs rvec position to cgal Point class */
    std::vector<Point_3> waterPoints = fromGmxtoCgalPosition<Point_3>(waterCoordinates, 3);
    
    waterDataHandle.startFrame(frnr, fr.time);
    volumeDataHandle.startFrame(frnr, fr.time);

    
    
    for (unsigned int i = 0; i < selectionListAlpha.size(); i++)
    {	
	gmx::ConstArrayRef<rvec> alphaCoordinates = selectionListAlpha.at(i).coordinates();

	/* Create CGAL Point_3 vector */
	std::vector<Point_3> alphaPoints = fromGmxtoCgalPosition<Point_3>(alphaCoordinates);
	
	/* Vector of buried water */
	std::vector<int> buriedWaterVector;

	/* Alpha shape computation */
	alphaShapeModulePtr_->build(alphaPoints, alphaValue_);    
	buriedWaterVector = alphaShapeModulePtr_->locate(waterPoints, true);	
	waterDataHandle.setPoint(i, buriedWaterVector.size());
	volumeDataHandle.setPoint(i, alphaShapeModulePtr_->volume());

	waterPresence_.push_back(std::vector<bool>(waterPoints.size(), false));
	for (auto &water : buriedWaterVector)
	{
	    waterPresence_.back().at(water) = true;
	}
	
	
	if (frnr == 10)
	{
	    std::ofstream oss;
	    std::string ofs;
	    oss.open("surface.off");
	    alphaShapeModulePtr_->writeOff(ofs);
	    oss << ofs;
	    oss.close();
	    
	    oss.open("water.pos");
	    oss << "Frame\n";
	    for (auto &water : buriedWaterVector)
	    {
		oss << waterPoints.at(water) << "\n";
	    }
	    oss.close();
	}
     
    }
    
    /* Store the output */
    
    waterDataHandle.finishFrame();
    volumeDataHandle.finishFrame();
}


void AlphaShape::finishAnalysis(int nframes)
{
    std::vector<int> lifetime(selectionWater_.posCount()/3);
    std::vector<int> cumul(selectionWater_.posCount()/3);
    std::vector<std::vector<bool> > waterPresenceT = transpose_2D(waterPresence_);
    for (auto &frame : waterPresence_)
    {
	for (unsigned int i = 0; i < frame.size(); ++i)
	{
	    if (frame.at(i))
	    {
		lifetime.at(i) += 1;
	    }
	}
    }
    
    int count = 0;
    int positive = 0;
    double stds = 0;
    for (auto &water : lifetime)
    {
	if (water >= nframes)
	{
	    std::cout << water << " ";
	    count += water;
	}
	if (water > 0)
	{
	    positive++;

	}
    }
    double average = 1.0*count/positive;
    for (auto &water : lifetime)
    {
	if (water > 0)
	{
	    stds += pow(water-average, 2.0);
	}
    }
    stds /= positive-1;
    std::cout << std::endl;
    std::cout << average << "+- " << pow(stds, 0.5) << std::endl;
    
    // std::ofstream oss;
    // oss.open("buried.ndx");
    // gmx::ConstArrayRef<int> waterindices = selectionWater_.atomIndices();
    // for (unsigned int i = 0; i < selectionListAlpha_.size(); ++i)
    // {
    //     oss << "[ Unit_" << i << " ]\n";
    // 	int per_line = 0;
    // 	int count = 0;
    // 	for (unsigned int j = 0; j < lifetime_.at(i).size(); ++j)
    // 	{
    // 	    if ( per_line == 4 )
    // 	    {
    // 		per_line = 0;
    // 		oss << "\n";
    // 	    } 
    // 	    if ( cumulLifetime_.at(i).at(j) > numFrameValue_ )
    // 	    {
    // 	        oss << waterindices.at(3*j+1) << " "
    // 		<< waterindices.at(3*j+2) << " "
    // 		<< waterindices.at(3*j+3) << " ";
    // 		++per_line;
    // 		++count;
    // 	    }
    // 	}
    // 	std::cout << " >> " << count << " water written\n";
    //     oss << "\n";
    // }
    // oss << "\n";
    // oss.close();
}


void AlphaShape::writeOutput()
{

}

int main(int argc, char *argv[])
{
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<AlphaShape>(argc, argv);
}
