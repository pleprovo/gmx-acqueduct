

#include "waternetwork.hpp"
#include "gromacs/pbcutil/pbc.h"

#include <iostream>
#include <fstream>
#include <limits>
#include <math.h>


struct custom_comparator {
    bool operator()(const std::pair<int, int>& a,
                    const std::pair<int, int>& b) const
    {
        return less_comparator(std::minmax(a.first, a.second),
                               std::minmax(b.first, b.second));
    }

    std::less<std::pair<int, int>> less_comparator;
};

double switch_function(double r, double r_on, double r_off)
{
    double sw = 0;

    if ( r <= r_off and r >= r_on )
    {
	sw = pow(pow(r,2)-pow(r_off,2),2)*(pow(r_off,2)+2*pow(r,2)-3*pow(r_on,2))/pow(pow(r_off,2)-pow(r_on,2),3);
    }
    else if ( r < r_on)
    {
	sw = 1.0;
    }
    return sw;
}

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

HydrogenBond computeHB(const rvec acceptor, const rvec donor, const rvec hydrogen)
{
    static float C = 3855; /* epsilon*sigma^6*sqrt(2/3) */
    static float D = 738; /* epsilon*sigma^4*sqrt(2/3) */
    static float theta_on = 0.25;
    static float theta_off = 0.03015369;
    static float r_on = 5.5;
    static float r_off = 6.5;
    double energy = 0.0;
    rvec AD, DH;
    HydrogenBond out;
    rvec_sub(hydrogen, donor, DH);
    rvec_sub(acceptor, donor, AD);
    
    double r = 10.0*norm(AD);
    double angle = cos_angle(DH, AD);
    if ( cos_angle(DH, AD) < -0.52 )
    {
	double r_switch = switch_function(r, r_on, r_off);
	double theta_switch = 1-switch_function(pow(angle, 2), pow(theta_off, 2),
						pow(theta_on, 2)); 
	energy = ((C/pow(r, 6.0))-(D/pow(r, 4.0)))*pow(cos(angle), 4.0)*r_switch*theta_switch;
    }
    else
    {
	energy = 0.0;
    }
    
    out.length = r;
    out.angle = angle;
    out.energy = energy;
    
    return out;
}

void checkHB(const gmx::AnalysisNeighborhoodPair &pair,
	     const gmx::Selection &waterSelection,
	     const gmx::ConstArrayRef<int> &mappedIds,
	     std::vector<HydrogenBond> &edgeVector)
{
    auto donor = waterSelection.position(mappedIds.at(pair.testIndex())).x();
    auto hydrogen1 = waterSelection.position(mappedIds.at(pair.testIndex())+1).x();
    auto hydrogen2 = waterSelection.position(mappedIds.at(pair.testIndex())+2).x();
    auto acceptor = waterSelection.position(mappedIds.at(pair.refIndex())).x();

    HydrogenBond out1, out2;
    out1 = computeHB(acceptor, donor, hydrogen1);
    out2 = computeHB(acceptor, donor, hydrogen2);
    edgeVector.push_back(out1);
    edgeVector.push_back(out2);
							  
}

WaterNetwork::WaterNetwork()
    : TrajectoryAnalysisModule("waternetwork", "Water network analysis tool"),
      cutoff_(0.65)
{
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
    options->addOption(gmx::SelectionOption("alpha").store(&calpha_).required()
		       .defaultSelectionText("Protein")
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
    selection.initOriginalIdsToGroup(top, INDEX_RES);
    gmx::ConstArrayRef<int> mappedIds = selection.mappedIds();
    unsigned int numGroup = mappedIds.back();
    unsigned int numAtoms = mappedIds.size();

    for (unsigned int i = 0; i < numGroup; i++)
    {
	std::cout << *top->atoms.resinfo[i].name << ", ";
    }
    std::cout << std::endl;
}


void WaterNetwork::initAnalysis(const gmx::TrajectoryAnalysisSettings &settings,
			   const gmx::TopologyInformation        &top)
{
    /* Init Selection */

    solvent_.initOriginalIdsToGroup(top.topology(), INDEX_RES);
    nb_water_ = solvent_.posCount()/3;

    // makeDonorAcceptorLists(calpha_, top.topology());
    
    /* Init neihborsearch cutoff value */
    nb1_.setCutoff(cutoff_);
    nb2_.setCutoff(cutoff_);
    nb1_.setMode(gmx::AnalysisNeighborhood::SearchMode::eSearchMode_Grid);
    nb2_.setMode(gmx::AnalysisNeighborhood::SearchMode::eSearchMode_Grid);
    
    /* Create ConstArray of indices for oxygen atoms */
    oxygenIndices_.reserve(nb_water_);
    for (int i = 0; i < nb_water_; i++)
    {
    	oxygenIndices_.push_back(3*i);
    }

    std::cout << "THE LAST WATER : " << oxygenIndices_.back() << std::endl;
    std::cout << "NUMBER OF WATER : " << oxygenIndices_.size() << std::endl;

    double on = 5.5;
    double off = 6.5;
    for (int r = 0; r < 50; r++)
    {
	std::cout << switch_function(5+r*0.05 , on, off) << ' ';
    }
    std::cout << std::endl;

    on = 0.25;
    off = 0.0301;
    for (int r = 0; r < 100; r++)
    {
	std::cout << 1-switch_function(0.26-r*0.001 , off, on) << ' ';
    }
    std::cout << std::endl;
    
    /* Set the number of column to store time dependent data */
    data_.setColumnCount(0, 4);

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
    const gmx::Selection           &sourcesel = pdata->parallelSelection(source_);
    const gmx::Selection           &calphasel = pdata->parallelSelection(calpha_);
    const gmx::Selection           &sinksel = pdata->parallelSelection(sink_);
    const gmx::Selection           &watersel = pdata->parallelSelection(solvent_);
    
    gmx::ConstArrayRef<int> oxygenArrayRef(oxygenIndices_);

    std::vector<Point_3> oxygenVec = fromGmxtoCgalPosition<Point_3>(watersel.coordinates(), 3);
    
    Delaunay DT(oxygenVec.begin(), oxygenVec.end());
    CGAL_assertion( DT.number_of_vertices() == oxygenVec.size() );   

    for(Delaunay::Finite_edges_iterator ei=DT.finite_edges_begin();ei!=DT.finite_edges_end(); ei++){
	
	std::cout << ei->first->vertex( (ei->second+1)%3)  << " ";
	std::cout << ei->first->vertex( (ei->second+2)%3)  << std::endl;
    }
    
    // std::cout << DT.number_of_edges() << " " << FDT.number_of_edges() << std::endl;
    
    /* Find all hydrogen bonds */
    // gmx::AnalysisNeighborhoodPair pair;
    // gmx::AnalysisNeighborhoodPositions sourcePos(sourcesel);
    // gmx::AnalysisNeighborhoodPositions calphaPos(calphasel);
    // gmx::AnalysisNeighborhoodPositions sinkPos(sinksel);    
    // gmx::AnalysisNeighborhoodPositions waterPos(watersel);
    
    // gmx::AnalysisNeighborhoodPositions oxygenPos = waterPos.indexed(oxygenArrayRef);
    
    // gmx::AnalysisNeighborhoodSearch searchAlpha = nb1_.initSearch(pbc, oxygenPos);
    // gmx::AnalysisNeighborhoodPairSearch pairSearchAlpha = searchAlpha.startPairSearch(calphaPos);
    // std::set<int> neighborSet;
    // std::vector<int> neighborVec;
    // while (pairSearchAlpha.findNextPair(&pair))
    // {
    //     neighborSet.insert(oxygenArrayRef.at(pair.refIndex()));
    // }
    // // std::cout << oxygenArrayRef.size() << " " << neighborSet.size() << std::endl;
    // std::copy(neighborSet.begin(), neighborSet.end(), std::back_inserter(neighborVec));
    // gmx::ConstArrayRef<int> neighborArrayRef(neighborVec);
  
    // gmx::AnalysisNeighborhoodPositions neighborPos = waterPos.indexed(neighborArrayRef);
    // gmx::AnalysisNeighborhoodSearch search1 = nb1_.initSearch(pbc, neighborPos);
    // gmx::AnalysisNeighborhoodPairSearch pairSearch = search1.startPairSearch(neighborPos);
	
    // gmx::AnalysisNeighborhoodSearch search2 = nb2_.initSearch(pbc, neighborPos);    
    // gmx::AnalysisNeighborhoodPairSearch pairSearchSource = search2.startPairSearch(sourcePos);
    // gmx::AnalysisNeighborhoodPairSearch pairSearchSink = search2.startPairSearch(sinkPos);

    std::set<int> sourceEdges;
    std::set<int> sinkEdges;
    std::vector<HydrogenBond> HBVector;
    // GraphModule g(nb_water_+2);
    
    // while (pairSearch.findNextPair(&pair))
    // {
    // 	if (pair.refIndex() != pair.testIndex())
    // 	{
	    
    // 	}
    // }
    
    // while (pairSearchSource.findNextPair(&pair))
    // {
    // 	sourceEdges.insert(pair.refIndex());
    // }

    // while (pairSearchSink.findNextPair(&pair))
    // {
    // 	sinkEdges.insert(pair.refIndex());
    // }


    
    dh.startFrame(frnr, fr.time);
    dh.setPoint(0, HBVector.size());
    dh.setPoint(1, sourceEdges.size());
    dh.setPoint(2, sinkEdges.size());
    dh.setPoint(3, 0/*flow*/);
    dh.finishFrame();
    // g.clear();
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

