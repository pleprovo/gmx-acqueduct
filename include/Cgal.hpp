
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>

#include <CGAL/property_map.h>
#include <boost/iterator/zip_iterator.hpp>
#include <utility>

struct edge {
    int i;
    int j;
    float length;
    float angle;
    float energy;
};

namespace cgal {

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3 Point_3;
    typedef boost::tuple<Point_3,int> Point_and_int;
    
    struct Node
    {
	int id;
	Point_3 h1;
	Point_3 h2;
	std::shared_ptr<std::vector<Point_3> > hydrogens;
    };
    
    typedef CGAL::Alpha_shape_vertex_base_3<K>                   Vb;
    typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vb>      Vbh;
    typedef CGAL::Triangulation_vertex_base_with_info_3<Node, K> Vbi;
    
    typedef CGAL::Alpha_shape_cell_base_3<K>                    Fb;
    typedef CGAL::Triangulation_data_structure_3<Vbh,Fb>        Tds;
    typedef CGAL::Triangulation_data_structure_3<Vbi>           Tdsi;

    typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Delaunay;
    typedef CGAL::Delaunay_triangulation_3<K, Tdsi, CGAL::Fast_location > DelaunayWithInfo;
    
    typedef CGAL::Alpha_shape_3<Delaunay>                       Alpha_shape_3;
    typedef CGAL::Search_traits_adapter<Point_and_int,
					CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
					CGAL::Search_traits_3<K> >
    Traits;
    typedef CGAL::Kd_tree<Traits>            Kd_tree_3;
    typedef CGAL::Fuzzy_sphere<Traits>       Fuzzy_sphere_3;
    /*
     * Compute the Volume of the Alpha Shape
     */
    float  getVolume(Alpha_shape_3 &s);

    /*
     * Find Points inside the Alpha Shape 
     */
    std::vector<int> searchPoints(Alpha_shape_3 &s, std::vector<Point_3> points);

    /*
     * Find Points in the vicinity of Kd Tree 
     */
    std::vector<int> searchPoints(Kd_tree_3 &s, std::vector<Point_3> points);
    
    /*
     * Compute Something over the edges of a given Triangulation
     */
    std::vector<edge> analyseEdges(DelaunayWithInfo &g);
    
}
