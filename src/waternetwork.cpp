#define _USE_MATH_DEFINES


#include "waternetwork.hpp"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"

#include <iostream>
#include <fstream>
#include <chrono>
#include <future>

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

std::ostream& operator << (std::ostream& os, const Site& site)
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
	      std::vector<Site> &sites)
{
    std::set<std::string> atomNames {
	"NE", /* 1H */
	    "NH1", /* 2H */
	    "NH2", /* 2H */
	    "ND1", /* 1H */
	    "ND2", /* 2H */
	    "OD1", /* 1H */
	    "OD2", /* 1H */
	    "NE1", /* 1H */
     	    "NE2", /* 1H */
	    "OE1", /* 0H */
	    "OE2", /* 0H */
	    "NZ",  /* 3H */
	    "OG",  /* 1H */
	    "OG1", /* 1H */
	    "OH",  /* 1H */
	    "O",   /* 0H */
	    "N",   /* 1H */
	    "SG",  /* 1H */
	    "OW", /* 2H */
	    "OCD"};

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
	    
	    Site site;
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


// std::vector<std::pair<std::pair<Site_ptr, Site_ptr>, hydrogenBond > >
// calc_edges (const cgal::DelaunayWithInfo::Finite_edges_iterator &begin_it,
// 	    const cgal::DelaunayWithInfo::Finite_edges_iterator &end_it,
// 	    const std::shared_ptr<cgal::computeEdgeEnergy> computator)
// {
//     std::vector<std::pair<std::pair<Site_ptr, Site_ptr>,
// 			  hydrogenBond > > segments;

//     for ( cgal::DelaunayWithInfo::Finite_edges_iterator it = begin_it;
// 	  it != end_it; it++ )
//     {
// 	auto v1 = it->first->vertex(cgal::DelaunayWithInfo::cw(it->second));
// 	auto v2 = it->first->vertex(cgal::DelaunayWithInfo::ccw(it->second));
// 	cgal::Info s1 = v1->info();
// 	cgal::Info s2 = v2->info();
// 	bool isedge = false;
// 	if ( s1.first.type != s2.first.type ) { isedge = true; }
// 	if ( s1.first.type == 2 && s2.first.type == 2) { isedge = true; }

// 	if (isedge)
// 	{
// 	    hydrogenBond hb = computator->compute(it);
// 	    if (hb.energy > 0.0)
// 	    {
// 		segments.push_back(std::make_pair(
// 				       std::make_pair(
// 					   std::make_shared<Site>(s1.first),
// 					   std::make_shared<Site>(s2.first)), hb));
// 	    }
// 	}
//     }
//     return segments;
// }

std::vector< HydrogenBond > make_edges (std::vector<cgal::Node>& nodes,
					unsigned int start, unsigned int stop)
{
    std::vector< HydrogenBond > results;

    for (unsigned int i = start; i < stop; ++i)
    {
	for ( unsigned int j = 1; j < i; ++j )
	{
	    bool isedge = false;

	    cgal::Vector_3 OO = *nodes.at(j).point - *nodes.at(i).point;
	    float distance = CGAL::sqrt(OO*OO);
	    OO /= distance;
	    distance *= 10.0;
    
	    if ( nodes.at(i).site->type != nodes.at(j).site->type ) { isedge = true; }
	    if ( nodes.at(i).site->type == WATER and nodes.at(j).site->type == WATER )
	    {
		isedge = true;
	    }
	    if ( distance > 6.5 ) { isedge = false; } 
	    if ( nodes.at(i).hydrogens.empty() and nodes.at(j).hydrogens.empty() )
	    {
		isedge = false;
	    } 

	    if ( isedge )
	    {
		std::vector<float> energies;
		std::vector<float> coss;
		
		if ( !nodes.at(i).hydrogens.empty() )
		{
		    for ( cgal::Point_ptr& hydrogen : nodes.at(i).hydrogens )
		    {
			cgal::Vector_3 DH = *hydrogen - *nodes.at(i).point;
			DH /= CGAL::sqrt(DH*DH);
			float proj = DH*OO;
			energies.push_back(cgal::computeEnergy(distance, proj));
			coss.push_back(proj);
		    }
		}
		
		if ( !nodes.at(j).hydrogens.empty() )
		{
		    for ( cgal::Point_ptr& hydrogen : nodes.at(j).hydrogens )
		    {
			cgal::Vector_3 DH = *hydrogen - *nodes.at(j).point;
			DH /= CGAL::sqrt(DH*DH);
			float proj = DH*(-1.0*OO);
			energies.push_back(cgal::computeEnergy(distance, proj));
			coss.push_back(proj);
		    }
		}
		
		// float cosmax = -0.1786;
		float cosmax = 0.80;
		int pos = -1;
		for ( unsigned int i = 0; i < coss.size(); ++i )
		{
		    if ( coss.at(i) > cosmax )
		    {
			cosmax = coss.at(i);
			pos = i;
		    }
		}

		if ( pos != -1 )
		{
		    results.push_back(HydrogenBond{
			    std::make_pair(nodes.at(i).site, nodes.at(j).site),
				FORWARD,
				distance,
				coss.at(pos),
				energies.at(pos)});
		}
	    }
	}
    }
    return results;
}

int make_edges_parallel(std::vector< cgal::Node >& nodes,
			std::vector< HydrogenBond >& edges)
{
    /* Make chunks */
    int n = std::thread::hardware_concurrency();
    int chunk_size = nodes.size() / n;

    /* Make Futures */
    std::vector< std::future< std::vector< HydrogenBond> > > futures;

    /* Launch Futures */
    for ( int k = 0; k < n-1; k++ )
    {
    	futures.push_back(std::async(std::launch::async,
    				     make_edges,
    				     std::ref(nodes),
    				     chunk_size*k,
    				     chunk_size*(k+1)));
    }
    
    edges = make_edges(nodes, chunk_size*(n-1), nodes.size());

    /* Await results */
    for ( std::future< std::vector<int> >& fut : futures )
    {
    	fut.wait();
    }

    /* Get result */
    for ( std::future< std::vector<int>>& fut : futures )
    {
	fut.wait();
    	std::vector< HydrogenBond  >chunk = fut.get();
    	edges.insert(edges.end(),
    			std::make_move_iterator(chunk.begin()),
    			std::make_move_iterator(chunk.end()));
    }
    
    return edges.size();
}


WaterNetwork::WaterNetwork() : alphavalue_(0.65),
			       lengthon_(5.5),
			       lengthoff_(6.5),
			       angleon_(0.25),
			       angleoff_(0.0301),
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
    
    options->addOption(gmx::SelectionOption("select").store(&solvent_).required()
		       .defaultSelectionText("Water")
		       .description(""));    
    
    options->addOption(gmx::SelectionOption("protein").store(&protein_).required()
		       .defaultSelectionText("Protein")
		       .description(""));
    options->addOption(gmx::SelectionOption("alpha").store(&alpha_).required()
		       .defaultSelectionText("Backbone")
		       .description(""));    

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
    /* Make Protein Sites */
    // protein_.initOriginalIdsToGroup(top_, INDEX_RES);
    if (protein_.isValid())
    {
	if (makeSites(protein_, top_, proteinSites_) > 0)
	{
	    std::clog << "Protein has : " << proteinSites_.size()
		      << " Sites of out " << protein_.posCount() << " atoms."
		      << std::endl;
	}
    }

    /* Make Solvent Sites */
    if (solvent_.isValid())
    {
    	if (makeSites(solvent_, top_, solventSites_) > 0)
    	{
    	    std::clog << "Solvent has : " << solventSites_.size()
		      << " Sites of out " << solvent_.posCount() << " atoms."
		      << std::endl;
	    for ( Site& site : solventSites_ )
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
    for ( Site& site : proteinSites_ )
    {
	outputStream_ << site;
    }
    for ( Site& site : solventSites_ )
    {
	outputStream_ << site;
    }
    
    outputStream_.close();
    
    std::clog << "Total number of sites : "
	      << proteinSites_.size()+solventSites_.size() << std::endl;
    
    /* Set the number of column to store time dependent data */
    filterData_.setColumnCount(0, 2); 
    graphData_.setColumnCount(0, 6);

    /* Init the average module  */ 
    avem_.reset(new gmx::AnalysisDataAverageModule());
    filterData_.addModule(avem_);

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
    
    const gmx::Selection &proteinsel = pdata->parallelSelection(protein_);
    const gmx::Selection &solventsel = pdata->parallelSelection(solvent_);
    const gmx::Selection &alphasel = pdata->parallelSelection(alpha_);

    auto start = std::chrono::high_resolution_clock::now();
    
    /* Convert positions from gromacs rvec to cgal point *///////////////////////////
    std::vector<cgal::Point_3> alphaPoints
	= fromGmxtoCgalPosition<cgal::Point_3>(alphasel.coordinates());
    std::vector<cgal::Point_3> proteinPoints
	= fromGmxtoCgalPosition<cgal::Point_3>(proteinsel.coordinates());
    std::vector<cgal::Point_3> solventPoints
	= fromGmxtoCgalPosition<cgal::Point_3>(solventsel.coordinates());
    std::vector<cgal::Point_3> oxygenPoints
	= fromGmxtoCgalPosition<cgal::Point_3>(solventsel.coordinates(), 3);
    
    auto stop0 = std::chrono::high_resolution_clock::now(); 
    
    /* Filter Solvent Sites *////////////////////////////////////////////////////////
    CGAL::spatial_sort(alphaPoints.begin(), alphaPoints.end());
    cgal::Alpha_shape_3 alphaShape(alphaPoints.begin(), alphaPoints.end());
    alphaShape.set_alpha(1.0);
    std::vector<int> filteredPoints;
    int num_filtered = cgal::filterPointsParallel(oxygenPoints, alphaShape, filteredPoints);

    auto stop1 = std::chrono::high_resolution_clock::now();
    
    /* Make Nodes *////////////////////////////////////////////////////////////////
    std::vector< cgal::Node > nodes;
    nodes.reserve(proteinSites_.size()+filteredPoints.size());
    
    for (unsigned int i = 0; i < proteinSites_.size(); i++)
    {
	std::vector<std::shared_ptr<cgal::Point_3>> hydrogens;  
	for (unsigned int i = 0; i < proteinSites_.at(i).nbHydrogen; i++)
	{
	    hydrogens.push_back(std::make_shared<cgal::Point_3>(
				    proteinPoints.at(proteinSites_.at(i).atmIndex + i+1)));
	}

	nodes.push_back(cgal::Node{std::make_shared<cgal::Point_3>
		    (proteinPoints.at(proteinSites_.at(i).atmIndex)),
		    std::make_shared<Site>(proteinSites_.at(i)),
		    hydrogens});
    }

    for (unsigned int i =0; i < filteredPoints.size(); i++)
    {	
	int index = filteredPoints.at(i);
	std::vector<std::shared_ptr<cgal::Point_3>> hydrogens;
	for (unsigned int i = 1; i <= solventSites_.at(index).nbHydrogen; i++)
	{
	    hydrogens.push_back(std::make_shared<cgal::Point_3>(
				    solventPoints.at(index*3 + i)));
	}

	nodes.push_back(cgal::Node{std::make_shared<cgal::Point_3>
		    (oxygenPoints.at(index)),
		    std::make_shared<Site>(solventSites_.at(index)),
		    hydrogens});
    }
    
    int num_nodes = nodes.size();

    auto stop2 = std::chrono::high_resolution_clock::now();
    
    /* Make Edges *///////////////////////////////////////////////////////////////////

    std::vector< HydrogenBond > segments;
    
    int num_edges = make_edges_parallel(nodes, segments);
    
    auto stop3 = std::chrono::high_resolution_clock::now();

    /* Write Graph *///////////////////////////////////////////////////////////////////
    
    outputStream_ << "#--- frame " << frnr << " ---\n";
    for (unsigned int i = 0; i < segments.size(); i++)
    {
    	outputStream_ << std::setw(6) << segments.at(i).sites.first->id
    		      << std::setw(6) << segments.at(i).sites.second->id
    		      << std::setw(8) << segments.at(i).energy
    		      << std::setw(8) << segments.at(i).length
    		      << std::setw(8) << segments.at(i).angle
    		      << "\n";
	
    }
    outputStream_.flush();

    auto stop4 = std::chrono::high_resolution_clock::now();

    auto total = std::chrono::duration_cast<std::chrono::milliseconds>(stop4 - start); 
    auto convert = std::chrono::duration_cast<std::chrono::milliseconds>(stop0 - start);
    auto filter = std::chrono::duration_cast<std::chrono::milliseconds>(stop1 - stop0);
    auto fill = std::chrono::duration_cast<std::chrono::milliseconds>(stop2 - stop1); 
    auto connect = std::chrono::duration_cast<std::chrono::milliseconds>(stop3 - stop2);
    auto write = std::chrono::duration_cast<std::chrono::milliseconds>(stop4 - stop3);
    
    filterData.startFrame(frnr, fr.time);
    filterData.setPoint(0, num_filtered);
    filterData.setPoint(1, cgal::getVolume(alphaShape));
    filterData.finishFrame();
    
    graphData.startFrame(frnr, fr.time);
    graphData.setPoint(0, convert.count());  
    graphData.setPoint(1, filter.count());
    graphData.setPoint(2, fill.count());
    graphData.setPoint(3, connect.count());
    graphData.setPoint(4, write.count());
    graphData.setPoint(5, 1.0*num_edges/num_nodes);
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

