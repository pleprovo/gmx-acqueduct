  

#include "alphashape.hpp"

#include <iostream>
#include <fstream>


// Convert position from Gmx to CGAL type
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

struct edge_comparator {
    bool operator()(const std::pair<int, int> &a,
                    const std::pair<int, int> &b) const {
	return less_comparator(std::minmax(a.first, a.second),
			       std::minmax(b.first, b.second));
    }

    std::less<std::pair<int, int> > less_comparator;
};


AlphaShape::AlphaShape()
    : TrajectoryAnalysisModule("AlphaShape", "Protein Surface analysis")
{
    registerAnalysisDataset(&data_, "avepop");

    alphaShapeModule_ = std::make_shared<AlphaShapeModule>();
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

    options->addOption(gmx::FileNameOption("o")
		       .filetype(gmx::eftPlot).outputFile()
		       .store(&fnPopulation_).defaultBasename("alpha-stats")
		       .description("Collection of analysis properties through time"));
    
    options->addOption(gmx::FileNameOption("os")
		       .filetype(gmx::eftUnknown).outputFile()
	               .store(&fnSurface_).defaultBasename("surface.off")
		       .description("Write the OFF file of the last alpha shape"));
    
    options->addOption(gmx::SelectionOption("reference")
		       .store(&alphasel_).required()
		       /*.defaultSelectionText("Calpha")*/
		       .description("Reference group to calculate alpha shape (default C-alphas)"));
    options->addOption(gmx::SelectionOption("select")
		       .store(&watersel_).required()
		       .defaultSelectionText("Water")
		       .description("Groups to calculate graph properties (default Water)"));

    settings->setFlag(gmx::TrajectoryAnalysisSettings::efRequireTop);
    settings->setFlag(gmx::TrajectoryAnalysisSettings::efUseTopX);
}


void AlphaShape::initAnalysis(const gmx::TrajectoryAnalysisSettings &settings,
			 const gmx::TopologyInformation        &top)
{
    /* Init Selection */
    watersel_.initOriginalIdsToGroup(top.topology(), INDEX_RES);

    /* Set the number of column to store time dependent data */
    data_.setColumnCount(0, 3);

    /* Init the average module  */ 
    avem_.reset(new gmx::AnalysisDataAverageModule());
    data_.addModule(avem_);

    /* Init the Plot module for the time dependent data */
    if (!fnPopulation_.empty()) {
	gmx::AnalysisDataPlotModulePointer plotm(
	    new gmx::AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnPopulation_);
        plotm->setTitle("Average distance");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Distance (nm)");
        data_.addModule(plotm);
    }
}


void AlphaShape::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
			   gmx::TrajectoryAnalysisModuleData *pdata)
{
    gmx::AnalysisDataHandle         dh     = pdata->dataHandle(data_);
    const gmx::Selection           &alphasel = pdata->parallelSelection(alphasel_);
    const gmx::Selection           &watersel = pdata->parallelSelection(watersel_);

    /* Get water Position and indices */
    gmx::ConstArrayRef<rvec> waterCoordinates = watersel.coordinates();
    gmx::ConstArrayRef<rvec> alphaCoordinates = alphasel.coordinates();

    /* Create CGAL Point_3 vector */
    std::vector<Point_3> alphaPoints = fromGmxtoCgalPosition<Point_3>(alphaCoordinates);
    std::vector<Point_3> waterPoints = fromGmxtoCgalPosition<Point_3>(waterCoordinates, 3);

    /* Vector of buried water */
    std::vector<int> buriedWaterVector;

    /* Alpha shape computation */
    alphaShapeModule_->build(alphaPoints, 3.0);    
    buriedWaterVector = alphaShapeModule_->locate(waterPoints, TRUE);

    if (frnr == 1) {
	std::string outputString;
	alphaShapeModule_->writeOff(outputString);
	// std::ofstream out(fnSurface_);
	// out << outputString;
	// out.close();
	
	std::ofstream OutFile;
	OutFile.open(fnSurface_, std::ios::out | std::ios::binary);
	OutFile.write( (char*)&outputString, sizeof(outputString));
	OutFile.close();
    }
    
    /* Store the output */
    dh.startFrame(frnr, fr.time);
    dh.setPoint(0, alphaShapeModule_->volume());
    dh.setPoint(1, buriedWaterVector.size());
    dh.setPoint(2, 0.0);
    dh.finishFrame();
    
}


void AlphaShape::finishAnalysis(int /*nframes*/)
{

}


void AlphaShape::writeOutput()
{

}

int main(int argc, char *argv[])
{
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<AlphaShape>(argc, argv);
}
