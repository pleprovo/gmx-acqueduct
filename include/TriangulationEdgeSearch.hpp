

#ifndef TRIANGULATIONEDGESEARCH_HPP
#define TRIANGULATIONEDGESEARCH_HPP

#include "Cgal.hpp"
#include <CGAL/Triangulation_hierarchy_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

struct Info
{
    int id;
    bool isSuperNode = false;
    std::shared_ptr<Point_3> h1;
    std::shared_ptr<Point_3> h2;
    std::shared_ptr<std::vector<Point_3> > hydrogens;
};

typedef CGAL::Triangulation_vertex_base_with_info_3<Info, K>    Vbi;
typedef CGAL::Triangulation_data_structure_3<Vbi>               Tdsi;
typedef CGAL::Delaunay_triangulation_3<K, Tdsi> DelaunayWithInfo;


class TriangulationEdgeSearch : public EdgeSearchInterface
{ 
public:
    std::vector<Edge> search(const std::vector<Node> nodes);
    
};

#endif
