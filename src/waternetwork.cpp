#define _USE_MATH_DEFINES


#include "waternetwork.hpp"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"

#include <iostream>
#include <fstream>
#include <chrono>
#include <future>

double switch_function_radius(const double r, const double r_on, const double r_off)
{
    double sw = 0.0;
		
    if ( r_off > r && r > r_on ) {
	sw = pow(pow(r,2)-pow(r_off,2),2);
	sw *= pow(r_off,2)+2*pow(r,2)-3*pow(r_on,2);
	sw /= pow(pow(r_off,2)-pow(r_on,2),3);	
    } else if ( r < r_on ) {
	sw = 1.0;
    }
    return sw;
}

double switch_function_angle(const double a, const double a_on, const double a_off)
{
    double sw = 0.0;
		
    if ( a_off < a && a < a_on ) {
	sw = pow(pow(a,2)-pow(a_off,2),2);
	sw *= pow(a_off,2)+2*pow(a,2)-3*pow(a_on,2);
	sw /= pow(pow(a_off,2)-pow(a_on,2),3);
	sw *= -1.0;
    } else if ( a > a_on ) {
	sw = 1.0;
    }
    return sw;
}

double computeEnergy(const double r,
		     const double cosine,
		     const double r_on = 5.5,
		     const double r_off = 6.5,
		     const double theta_on = 0.25,
		     const double theta_off = 0.0)
{
    static float C = 3855; /* epsilon*sigma^6*sqrt(2/3) */
    static float D = 738; /* epsilon*sigma^4*sqrt(2/3) */
    double radius_switch = switch_function_radius(r*r, r_on*r_on, r_off*r_off);
    double angle_switch = switch_function_angle(cosine*cosine,
						theta_off, theta_on); 
    double energy = ((C/pow(r, 6.0))-(D/pow(r, 4.0)));
    energy *= pow(cosine, 4.0);
    energy *= radius_switch;
    energy *= angle_switch;
    return energy;
}

template <typename T>
std::vector<T> fromGmxtoCgalPosition(const gmx::ArrayRef<const rvec>& coordinates,
				     const int increment=1)
{
    std::vector<T> cgalPositionVector;
    cgalPositionVector.reserve(coordinates.size());
    for (unsigned int i = 0; i < coordinates.size(); i += increment)
    {
	cgalPositionVector.push_back(T(coordinates.at(i)[XX],
				       coordinates.at(i)[YY],
				       coordinates.at(i)[ZZ]));
    }
    return cgalPositionVector;
}

std::ostream& operator << (std::ostream& os, const SiteInfo& site)
{
    os << std::setw(7) << site.id 
       << std::setw(7) << site.atmIndex 
       << std::setw(7) << site.name 
       << std::setw(7) << site.resIndex 
       << std::setw(7) << site.resName 
       << std::setw(7) << site.type 
       << std::setw(7) << site.nbHydrogen << "\n";
    return os;
}

int makeSites(gmx::Selection &selection,
	      const gmx::TopologyInformation *top,
	      std::vector<SiteInfo> &sites)
{
    std::set<std::string> atomNames {
	"OW",  /* WAT 2H */	
	    "O",   /* BCK 0H */
	    "N",   /* BCK 1H */
	    "SG",  /* CYS 1H */
     	    "NE2", /* GLN 2H */
	    "OE1", /* GLN 0H */
	    "OD1", /* ASN 0H */
	    "ND2", /* ASN 2H */
	    "OG",  /* SER 1H */
	    "OG1", /* THR 1H */
	    "OH",  /* TYR 1H */
	    "NE1", /* TRP 1H */
	    "NZ",  /* LYS 3H */
	    "NH1", /* ARG 2H */
	    "NH2", /* ARG 2H */
	    "ND1", /* HIS 1H */
	    "NE2", /* HIS 1H */};
    
    std::map<std::string, SiteType> siteType;
    siteType["NE"] = DONOR; /* NE */
    siteType["NH1"] = DONOR; /* NH1 */
    siteType["NH2"] = DONOR; /* NH2 */
    siteType["ND1"] = DONOR; /* ND1 */
    siteType["ND2"] = DONOR; /* ND2 */
    siteType["OD1"] = DONOR; /* OD1 */
    siteType["OD2"] = DONOR; /* OD2 */
    siteType["NE1"] = DONOR; /* NE1 */
    siteType["NE2"] = DONOR; /* NE2 */
    siteType["OE1"] = ACCEPTOR; /* OE1 */
    siteType["OE2"] = ACCEPTOR; /* OE2 */
    siteType["NZ"] = DONOR; /* NZ */
    siteType["OG"] = DONOR; /* OG */
    siteType["OH"] = DONOR; /* OH */
    siteType["O"] = ACCEPTOR; /* O */
    siteType["N"] = DONOR; /* N */
    siteType["SG"] = DONOR; /* SG */ 
    siteType["OW"] = WATER; /* OW */
    siteType["OCD"] = ACCEPTOR;

    std::string watres("WAT");
    
    gmx::ArrayRef<const int> indices = selection.atomIndices();

    unsigned int count = 0;
    for (unsigned int i = 0; i < indices.size(); i++)
    {
	int index = indices.at(i);
	
    	std::string name(*top->atoms()->atomname[index]);
	if ( atomNames.find(name) != atomNames.end() )
	{
	    count++;	    
	    unsigned int nhb = 0;
	    while (true)
	    {
	    	if ( index + nhb + 1 > indices.size() ) { break; }
	    	std::string probe(top->atoms()->atom[index + nhb + 1].elem);
	    	if ( probe.find("H") == std::string::npos ) { break; }
	    	++nhb;
	    }
	    
	    SiteInfo site;
	    site.id = count;
	    site.name = name;
	    site.atmIndex = index;
	    site.resIndex = top->atoms()->atom[index].resind;
	    site.nbHydrogen = nhb;
	    site.type = siteType[site.name];
	    if (watres.compare(*top->atoms()->resinfo[site.resIndex].name) == 0)
	    {
	    	site.type = WATER;
	    }
	    site.resName = std::string(*top->atoms()->resinfo[site.resIndex].name);
	    sites.push_back(site);
	}
    }

    return sites.size();
}

std::vector<HydrogenBond> make_edges(const std::vector<std::pair<int, int>>::iterator& begin,
				     const std::vector<std::pair<int, int>>::iterator& end,
				     const std::vector<Point>& sitePoints,
				     const std::vector<std::vector<Point_ptr>>& hydrogens,
				     const std::vector<SiteInfo_ptr>& sites)
{
    std::vector< HydrogenBond > results;

    for (std::vector<std::pair<int, int>>::iterator it = begin; it != end; ++it)
    {
	auto pair = *it;
	bool isedge = false;
	Vector_3 OO = sitePoints.at(pair.first) - sitePoints.at(pair.second);
	float distance = CGAL::sqrt(OO*OO);
	OO /= distance;
	distance *= 10.0;
    
	if ( sites.at(pair.first)->type != sites.at(pair.second)->type ) { isedge = true; }
	if ( sites.at(pair.first)->type == WATER and sites.at(pair.second)->type == WATER )
	{
	    isedge = true;
	}
	if ( distance > 6.5 ) { isedge = false; } 
	if ( hydrogens.at(pair.first).empty() and hydrogens.at(pair.second).empty() )
	{
	    isedge = false;
	} 

	if ( isedge )
	{
	    std::vector<float> energies;
	    std::vector<float> coss;
		
	    if ( !hydrogens.at(pair.first).empty() )
	    {
		for ( const Point_ptr& hydrogen : hydrogens.at(pair.first) )
		{
		    Vector_3 DH = *hydrogen - sitePoints.at(pair.first);
		    DH /= CGAL::sqrt(DH*DH);
		    float proj = DH*OO;
		    energies.push_back(computeEnergy(distance, proj));
		    coss.push_back(proj);
		}
	    }
		
	    if ( !hydrogens.at(pair.second).empty() )
	    {
		for ( const Point_ptr& hydrogen : hydrogens.at(pair.second) )
		{
		    Vector_3 DH = *hydrogen - sitePoints.at(pair.second);
		    DH /= CGAL::sqrt(DH*DH);
		    float proj = DH*(-1.0*OO);
		    energies.push_back(computeEnergy(distance, proj));
		    coss.push_back(proj);
		}
	    }
		
	    // float cosmax = -0.1786;
	    float cosmax = 0.0;
	    int pos = -1;
	    for ( unsigned int i = 0; i < coss.size(); ++i )
	    {
		if ( coss.at(i) > cosmax )
		{
		    cosmax = coss.at(i);
		    pos = i;
		}
	    }

	    if ( pos > 0 )
	    {
		results.push_back(HydrogenBond{
			SiteInfoPair(sites.at(pair.first), sites.at(pair.second)),
			    FORWARD,
			    distance,
			    coss.at(pos),
			    energies.at(pos)});
	    }
	}
    }
    return results;
}

// int make_edges_parallel(const std::vector<std::pair<int, int>>& pairs,
// 			const std::vector<Point>& sitePoints,
// 			const std::vector<std::vector<Point_ptr>>& hydrogens,
// 			const std::vector<SiteInfo_ptr>& sites)
// {
//     /* Make chunks */
//     int n = std::thread::hardware_concurrency();
//     int chunk_size = pairs.size() / n;

//     /* Make Futures */
//     std::vector< std::future< std::vector< HydrogenBond> > > futures;

//     /* Launch Futures */
//     for ( int k = 0; k < n-1; k++ )
//     {
//     	futures.push_back(std::async(std::launch::async,
//     				     make_edges,
//     				     std::ref(nodes),
//     				     chunk_size*k,
//     				     chunk_size*(k+1)));
//     }
    
//     edges = make_edges(nodes, chunk_size*(n-1), nodes.size());

//     /* Await results */
//     for ( std::future< std::vector<int> >& fut : futures )
//     {
//     	fut.wait();
//     }

//     /* Get result */
//     for ( std::future< std::vector<int>>& fut : futures )
//     {
// 	fut.wait();
//     	std::vector< HydrogenBond  >chunk = fut.get();
//     	edges.insert(edges.end(),
//     			std::make_move_iterator(chunk.begin()),
//     			std::make_move_iterator(chunk.end()));
//     }
    
//     return edges.size();
// }


WaterNetwork::WaterNetwork() : alphaValue_(0.65),
			       lengthOn_(5.5),
			       lengthOff_(6.5),
			       angleOn_(0.25),
			       angleOff_(0.0301),
			       top_(nullptr)
{
    registerAnalysisDataset(&filterData_, "filter");
    registerAnalysisDataset(&graphData_, "graph");
}


void WaterNetwork::initOptions(gmx::IOptionsContainer          *options,
			       gmx::TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "This is a template for writing your own analysis tools for",
	"GROMACS. The advantage of using GROMACS for this is that you",
        "have access to all information in the topology, and your",
        "program will be able to handle all types of coordinates and",
        "trajectory files supported by GROMACS. In addition,",
        "you get a lot of functionality for free from the trajectory",
        "analysis library, including support for flexible dynamic",        "selections. Go ahead an try it![PAR]",
        "To get started with implementing your own analysis program,",
        "follow the instructions in the README file provided.",
        "This template implements a simple analysis programs that calculates",
        "average distances from a reference group to one or more",
        "analysis groups."
    };

    settings->setHelpText(desc);

    options->addOption(gmx::FileNameOption("o").filetype(gmx::eftPlot).outputFile()
		       .store(&fnFilter_).defaultBasename("filter")
		       .description("Collection of analysis properties through time"));

    options->addOption(gmx::FileNameOption("og").filetype(gmx::eftPlot).outputFile()
		       .store(&fnGraph_).defaultBasename("graph")
		       .description("Collection of analysis properties through time"));
    
    options->addOption(gmx::SelectionOption("select").store(&solventSel_).required()
		       .defaultSelectionText("Water")
		       .description(""));    
    
    options->addOption(gmx::SelectionOption("protein").store(&proteinSel_).required()
		       .defaultSelectionText("Protein")
		       .description(""));
    options->addOption(gmx::SelectionOption("alpha").store(&alphaSel_).required()
		       .defaultSelectionText("Backbone")
		       .description(""));    

    options->addOption(gmx::DoubleOption("alphavalue").store(&alphaValue_)
		       .description("Alpha value for alpha shape computation (default 1.0nm)"));
    options->addOption(gmx::DoubleOption("don").store(&lengthOn_)
		       .description("Distance Cutoff for energy calculation (default 0.25 nm)"));
    options->addOption(gmx::DoubleOption("doff").store(&lengthOff_)
		       .description("Distance Cutoff for energy calculation (default 0.25 nm)"));
    options->addOption(gmx::DoubleOption("aon").store(&angleOn_)
		       .description("Angle cutoff for energy calculation (default cos(90))"));
    options->addOption(gmx::DoubleOption("aoff").store(&angleOff_)
		       .description("Angle cutoff for energy calculation (default cos(90))")); 

}

void WaterNetwork::optionsFinished(gmx::TrajectoryAnalysisSettings *settings)
{

    settings->setFlag(gmx::TrajectoryAnalysisSettings::efRequireTop);
    settings->setFlag(gmx::TrajectoryAnalysisSettings::efUseTopX);

}

void WaterNetwork::initAnalysis(const gmx::TrajectoryAnalysisSettings &settings,
				const gmx::TopologyInformation        &top)
{
    top_ = &top;

    as_ = std::make_shared<AlphaShapeSurface>();
    as_->setAlphaValue(alphaValue_);

    mp_ = std::make_shared<DelaunayFindPairs>();
    
    /* Make Protein Sites */
    if (proteinSel_.isValid())
    {
	if (makeSites(proteinSel_, top_, proteinSites_) > 0)
	{
	    std::clog << "Protein has : " << proteinSites_.size()
		      << " Sites of out " << proteinSel_.posCount() << " atoms."
		      << std::endl;
	}
    }

    /* Make Solvent Sites */
    if (solventSel_.isValid())
    {
    	if (makeSites(solventSel_, top_, solventSites_) > 0)
    	{
    	    std::clog << "Solvent has : " << solventSites_.size()
		      << " Sites of out " << solventSel_.posCount() << " atoms."
		      << std::endl;
	    for ( SiteInfo& site : solventSites_ )
	    {
		site.id += proteinSites_.size()+1;
	    }
    	}
    }   

    
    
    outputStream_.open("node_list.txt");
    outputStream_ << std::setw(7) << "#Id" 
		  << std::setw(7) << "AtmId" 
		  << std::setw(7) << "Name" 
		  << std::setw(7) << "ResId" 
		  << std::setw(7) << "Res" 
		  << std::setw(7) << "Type" 
		  << std::setw(7) << "Hyd\n";
    outputStream_.flush();
    for ( SiteInfo& site : proteinSites_ )
    {
	outputStream_ << site;
    }
    for ( SiteInfo& site : solventSites_ )
    {
	outputStream_ << site;
    }
    
    outputStream_.close();
    
    std::clog << "Total number of sites : "
	      << proteinSites_.size()+solventSites_.size() << std::endl;
    
    /* Set the number of column to store time dependent data */
    filterData_.setColumnCount(0, 2); 
    graphData_.setColumnCount(0, 6);

    /* Init the Plot module for the time dependent data */
    /* Solvent filtering data */
    if (!fnFilter_.empty())
    {
	gmx::AnalysisDataPlotModulePointer plotm(
	    new gmx::AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnFilter_);
        plotm->setTitle("Filter Statistics");
        plotm->setXAxisIsTime();
        plotm->setYLabel("#");
        filterData_.addModule(plotm);
    }

    /* Init the Plot module for the time dependent data */
    /* Network Data */
    if (!fnGraph_.empty())
    {
	gmx::AnalysisDataPlotModulePointer plotm(
	    new gmx::AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnGraph_);
        plotm->setTitle("Graph Statistics");
        plotm->setXAxisIsTime();
        plotm->setYLabel("#");
        graphData_.addModule(plotm);
    }

    // Edge list file open
    outputStream_.open("edge_list.txt");
    outputStream_ << std::fixed << std::setprecision(4);
    outputStream_ << std::setw(6) << "#Id"
		  << std::setw(8) << "AtmId"
		  << std::setw(6) << "ResId"
		  << std::setw(6) << "Id"
		  << std::setw(8) << "AtmId"
		  << std::setw(6) << "ResId"
		  << std::setw(8) << "Energy"
		  << std::setw(8) << "Length"
		  << std::setw(8) << "Cos\n";
    outputStream_.flush();
}


void WaterNetwork::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
				gmx::TrajectoryAnalysisModuleData *pdata)
{
    gmx::AnalysisDataHandle filterData = pdata->dataHandle(filterData_);
    gmx::AnalysisDataHandle graphData = pdata->dataHandle(graphData_);
    
    const gmx::Selection &proteinsel = pdata->parallelSelection(proteinSel_);
    const gmx::Selection &solventsel = pdata->parallelSelection(solventSel_);
    const gmx::Selection &alphasel = pdata->parallelSelection(alphaSel_);

    auto start = std::chrono::high_resolution_clock::now();
    
    /* Convert positions from gromacs rvec to cgal point *///////////////////////////
    std::vector<Point> alphaPoints
	= fromGmxtoCgalPosition<Point>(alphasel.coordinates());
    std::vector<Point> proteinPoints
	= fromGmxtoCgalPosition<Point>(proteinsel.coordinates());
    std::vector<Point> solventPoints
	= fromGmxtoCgalPosition<Point>(solventsel.coordinates());
    std::vector<Point> oxygenPoints
	= fromGmxtoCgalPosition<Point>(solventsel.coordinates(), 3);
    
    auto stop0 = std::chrono::high_resolution_clock::now(); 
    
    /* Filter Solvent Sites *////////////////////////////////////////////////////////
    std::vector<int> filteredPoints;
    as_->make(alphaPoints);
    
    auto stop1 = std::chrono::high_resolution_clock::now();
    
    int num_filtered = as_->locate(oxygenPoints, filteredPoints);
    
    auto stop2 = std::chrono::high_resolution_clock::now();
    
    /* Make Nodes *////////////////////////////////////////////////////////////////
    std::vector<Point> sitePoints;
    std::vector<std::vector<Point_ptr>> hydrogenPoints;
    std::vector<SiteInfo_ptr> siteInfos;
    
    sitePoints.reserve(proteinSites_.size() + filteredPoints.size());
    hydrogenPoints.reserve(proteinSites_.size() + filteredPoints.size());
    siteInfos.reserve(proteinSites_.size() + filteredPoints.size());
    
    
    for (SiteInfo& info : proteinSites_)
    {
    	std::vector<Point_ptr> hydrogens(info.nbHydrogen);
	if (info.nbHydrogen != 0)
	{
	    for (unsigned int i = 0; i < info.nbHydrogen; i++)
	    {
		hydrogens.at(i) = std::make_shared<Point>(proteinPoints.at(info.id+i+1));
	    }
	}
	hydrogenPoints.push_back(hydrogens);
	siteInfos.push_back(std::make_shared<SiteInfo>(info));
    	sitePoints.push_back(proteinPoints.at(info.atmIndex));
    }

    for (int index : filteredPoints)
    {	
    	std::vector<Point_ptr> hydrogens(solventSites_.at(index).nbHydrogen);
    	if (solventSites_.at(index).nbHydrogen)
    	{
    	    for (unsigned int i = 0; i < solventSites_.at(index).nbHydrogen; i++)
    	    {
    		hydrogens.at(i) = std::make_shared<Point>(solventPoints.at(3*index+i));
    	    }
    	}
    	hydrogenPoints.push_back(hydrogens);
    	siteInfos.push_back(std::make_shared<SiteInfo>(solventSites_.at(index)));
    	sitePoints.push_back(solventPoints.at(3*index));
    }
    
    int num_sites = sitePoints.size();

    auto stop3 = std::chrono::high_resolution_clock::now();
    
    // /* Find Edges *///////////////////////////////////////////////////////////////////

    std::vector<std::pair<int, int>> pairs = mp_->find(sitePoints);
    
    auto stop4 = std::chrono::high_resolution_clock::now();

    // /* Find Hydrogen Bonds *//////////////////////////////////////////////////////////
    //TODO Make this function async

    /* Make chunks */
    int n = std::thread::hardware_concurrency();
    int chunksize = pairs.size() / n;
    // std::cout << " Concurrency " << n << " chunksize " << chunksize << std::endl;
    /* Make Futures */
    std::vector<std::future<std::vector<HydrogenBond>>> futures;

    /* Launch Futures */
    for ( int k = 0; k < n-1; k++ )
    {
    	futures.push_back(std::async(std::launch::async,
    				     std::ref(make_edges),
    				     pairs.begin()+chunksize*k,
    				     pairs.begin()+chunksize*(k+1),
    				     std::ref(sitePoints),
    				     std::ref(hydrogenPoints),
    				     std::ref(siteInfos)));
    }
    
    futures.push_back(std::async(std::launch::async,
				 std::ref(make_edges),
				 pairs.begin()+chunksize*(n-1),
				 pairs.end(),
				 std::ref(sitePoints),
				 std::ref(hydrogenPoints),
				 std::ref(siteInfos)));

    /* Get result */
    std::vector<std::vector<HydrogenBond>> results;
    for ( auto &fut : futures )
    {
	results.push_back(fut.get());
    	// std::vector<HydrogenBond> chunk = fut.get();
    }
    
    std::vector<HydrogenBond> edges;
    for ( auto &chunk : results )
    {
    	edges.insert(edges.end(),
    		     std::make_move_iterator(chunk.begin()),
    		     std::make_move_iterator(chunk.end()));
    }
    
    // std::vector<HydrogenBond> edges = make_edges(pairs.begin(), pairs.end(), sitePoints,
    // 						 hydrogenPoints, siteInfos);
    
    int num_edges = edges.size();

    auto stop5= std::chrono::high_resolution_clock::now();

    // /* Write Graph *///////////////////////////////////////////////////////////////////
    
    outputStream_ << "#--- frame " << frnr << " ---\n";
    for (unsigned int i = 0; i < edges.size(); i++)
    {
    	outputStream_ << std::setw(6) << edges.at(i).sites.first->id
    		      << std::setw(6) << edges.at(i).sites.second->id
    		      << std::setw(8) << edges.at(i).energy
    		      << std::setw(8) << edges.at(i).length
    		      << std::setw(8) << edges.at(i).angle
    		      << "\n";
	
    }
    outputStream_.flush();

    auto total = std::chrono::duration_cast<std::chrono::milliseconds>(stop4 - start); 
    auto convert = std::chrono::duration_cast<std::chrono::milliseconds>(stop0 - start);
    auto alpha = std::chrono::duration_cast<std::chrono::milliseconds>(stop1 - stop0);
    auto locate = std::chrono::duration_cast<std::chrono::milliseconds>(stop2 - stop1); 
    auto connect = std::chrono::duration_cast<std::chrono::milliseconds>(stop3 - stop2);
    auto energy = std::chrono::duration_cast<std::chrono::milliseconds>(stop4 - stop3);
    auto write = std::chrono::duration_cast<std::chrono::milliseconds>(stop5 - stop4);
    
    filterData.startFrame(frnr, fr.time);
    filterData.setPoint(0, num_filtered);
    filterData.setPoint(1, 0.0/*cgal::getVolume(alphaShape)*/);
    filterData.finishFrame();
    
    graphData.startFrame(frnr, fr.time);
    graphData.setPoint(0, 1.0*convert.count()/total.count()/*num_sites*/);  
    graphData.setPoint(1, 1.0*alpha.count()/total.count()/*num_edges*/);
    graphData.setPoint(2, 1.0*locate.count()/total.count()/*1.0*num_edges/num_sites*/);
    graphData.setPoint(3, 1.0*connect.count()/total.count()/*0.0*/);
    graphData.setPoint(4, 1.0*energy.count()/total.count()/*0.0*/);
    graphData.setPoint(5, 1.0*write.count()/total.count()/*0.0*/);
    graphData.finishFrame();

    // std::cout << "Num_nodes = " << num_nodes << " Num_edges = " << num_edges << std::endl;
}


void WaterNetwork::finishAnalysis(int /*nframes*/)
{
    outputStream_.close();
}


void WaterNetwork::writeOutput()
{
    
}


int main(int argc, char *argv[])
{
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<WaterNetwork>(argc, argv);
}

