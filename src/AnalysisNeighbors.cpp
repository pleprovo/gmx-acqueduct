
#include "Cgal.hpp"
#include "AnalysisNeighbors.hpp"


void AnalysisNeighbors::initialize(Config config)
{
    
}

Results AnalysisNeighbors::execute(const Frame &frame)
{

    /*
     *Converstion of positions set to cgal point vectors
     */
    std::vector<cgal::Point_3> points = fromGmxtoCgalPosition<cgal::Point_3>(frame.points.coordinates());
    std::vector<cgal::Point_3> proteins = fromGmxtoCgalPosition<cgal::Point_3>(frame.protein.coordinates());
    std::vector<cgal::Point_3> solvent = fromGmxtoCgalPosition<cgal::Point_3>(frame.solvent.coordinates());
    std::vector<cgal::Point_3> oxygens = fromGmxtoCgalPosition<cgal::Point_3>(frame.solvent.coordinates(), 3);

    std::vector<int> indices(oxygens.size());
    std::iota (std::begin(indices), std::end(indices), 0);

    cgal::Kd_tree_3 tree(
	boost::make_zip_iterator(boost::make_tuple(oxygens.begin(), indices.begin() )),
	boost::make_zip_iterator(boost::make_tuple(oxygens.end(), indices.end() ) )  
	);
    // cgal::Kd_tree_3 tree(oxygens.begin(), oxygens.end());
    
    std::vector<int> selected = cgal::searchPoints(tree, points);
    
    cgal::DelaunayWithInfo DTI;
    for (auto id : selected)
    {
    	cgal::DelaunayWithInfo::Vertex_handle vertex =  DTI.insert(oxygens.at(id));
    	vertex->info() = cgal::Node{3*id, solvent.at(3*id+1), solvent.at(3*id+2)}; 	
    }
    
    std::vector<edge> edges = cgal::analyseEdges(DTI);
    std::cout << DTI << std::endl;
    Results results{DTI.number_of_edges(), edges.size()};     
    return results;
}
