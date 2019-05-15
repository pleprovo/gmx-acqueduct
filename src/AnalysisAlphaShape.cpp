
#include "Cgal.hpp"
#include "AnalysisAlphaShape.hpp"

#include <map>

void StrategyAlpha::initialize(Config config)
{
    
}

Results StrategyAlpha::execute(const Frame &frame)
{

    /*
     *Converstion of positions set to cgal point vectors
     */
    std::vector<cgal::Point_3> points = fromGmxtoCgalPosition<cgal::Point_3>(frame.points.coordinates());
    std::vector<cgal::Point_3> proteins = fromGmxtoCgalPosition<cgal::Point_3>(frame.protein.coordinates());
    std::vector<cgal::Point_3> solvent = fromGmxtoCgalPosition<cgal::Point_3>(frame.solvent.coordinates());
    std::vector<cgal::Point_3> oxygens = fromGmxtoCgalPosition<cgal::Point_3>(frame.solvent.coordinates(), 3);

    // make maps between point_3 and sites
    
    std::map<cgal::Point_3, Site> siteMap;
    for (auto site : frame.solventSites)
    {
    	siteMap[solvent.at(site.index)] = site;
    }
    // std::map<cgal::Point_3, Site> proteinMap;
    
    
    cgal::Alpha_shape_3 alphaShape(points.begin(), points.end());
    alphaShape.set_alpha(1.0);
    std::vector<int> selected = cgal::searchPoints(alphaShape, oxygens);
    
    cgal::DelaunayWithInfo DTI;
    // for (auto site : frame.proteinSites)
    // {
    // 	siteMap[proteins.at(site.index)] = site;
    // 	cgal::DelaunayWithInfo::Vertex_handle vertex = DTI.insert(proteins.at(site.index));
    // 	vertex->info() = cgal::Node{site.index};
    // }
    // for (auto id : selected)
    // {
    // 	cgal::DelaunayWithInfo::Vertex_handle vertex =  DTI.insert(oxygens.at(id));
    // 	vertex->info() = cgal::Node{3*id, solvent.at(3*id+1), solvent.at(3*id+2)}; 	
    // }

    for (unsigned int i = 0; i < oxygens.size(); i++)
    {
    	cgal::DelaunayWithInfo::Vertex_handle vertex =  DTI.insert(oxygens.at(i));
    	vertex->info() = cgal::Node{3*i, solvent.at(3*i+1), solvent.at(3*i+2)};
    }
    
    std::vector<edge> edges = cgal::analyseEdges(DTI);
    Results results{selected.size(), edges.size(), cgal::getVolume(alphaShape)};
    return results;
    
}
