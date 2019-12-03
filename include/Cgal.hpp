
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/property_map.h>

// #include <boost/iterator/zip_iterator.hpp>

#include <tuple>

#include "utils.hpp"

enum NodeIndex {POINT, SITE, HYDROGENS};

namespace cgal {

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3 Point_3;
    typedef K::Vector_3 Vector_3;
    
    using Point_ptr = std::shared_ptr<Point_3>;
    using Info = std::pair<Site_ptr, std::vector<Point_ptr> >;

    struct Node
    {
	Point_ptr point;
	Site_ptr site;
	std::vector<Point_ptr> hydrogens;
    };

    using NodeTuple = std::tuple<Point_ptr, Site_ptr, std::vector<Point_ptr> >;
    
    typedef boost::tuple<Point_3, Info> Point_and_Info;
    
    typedef CGAL::Alpha_shape_vertex_base_3<K> Vb;

    // Definition for Alpha Shape 3
    typedef CGAL::Triangulation_data_structure_3<
	CGAL::Triangulation_hierarchy_vertex_base_3<Vb>,
	CGAL::Alpha_shape_cell_base_3<K>,
	CGAL::Parallel_tag> Tds;
    typedef CGAL::Delaunay_triangulation_3<K,
					   Tds,
					   CGAL::Fast_location> Delaunay;
    typedef CGAL::Alpha_shape_3<Delaunay> Alpha_shape_3;

    // Definition for Delaunay with info
    typedef CGAL::Triangulation_data_structure_3<
	CGAL::Triangulation_vertex_base_with_info_3<Info, K>,
	CGAL::Delaunay_triangulation_cell_base_3<K>,
	CGAL::Parallel_tag> Tdsi;
    typedef CGAL::Delaunay_triangulation_3<K,
					   Tdsi,
					   CGAL::Fast_location > DelaunayWithInfo;
    
    // Definition for Kd_tree
    typedef CGAL::Search_traits_adapter<
	Point_and_Info,
	CGAL::Nth_of_tuple_property_map<0, Point_and_Info>,
	CGAL::Search_traits_3<K> > Traits;
    
    typedef CGAL::Kd_tree<Traits>            Kd_tree_3;
    typedef CGAL::Fuzzy_sphere<Traits>       Fuzzy_sphere_3;


    /*
     * Distance Switch function 
     */
    double switch_function_radius(const double r, const double r_on, const double r_off);

    /*
     * Angle Switch function 
     */
    double switch_function_angle(const double a, const double a_on, const double a_off);

    /*
     * Hydrogen bond energy function
     */

    double computeEnergy(const double r,
			 const double cosine,
			 const double r_on = 5.5,
			 const double r_off = 6.5,
			 const double theta_on = 0.25,
			 const double theta_off = 0.0301);

    /*
     * Compute the volume of and Alpha Shape
     */
    float getVolume(Alpha_shape_3 &s);

    /* 
     * Find point in Alpha Shape 
     */
    std::vector<int> filterPoints (std::vector<Point_3>& points,
				   Alpha_shape_3& s,
				   int start,
				   int stop);

    /* 
     * Asynchronous version of filterPoints 
     */
    int filterPointsParallel(std::vector<Point_3>& points,
			     Alpha_shape_3& s,
			     std::vector<int>& results);
 
    // std::vector<int> filterPoints(std::vector<Point_3> &points,
    // 				  Kd_tree_3 &s)
    // {
	
    // 	std::set<int> indiceLocated;
    // 	std::vector<cgal::Point_and_int> neighbors;
    // 	for (unsigned  int i = 0; i < points.size(); i++)
    // 	{
    // 	    Fuzzy_sphere_3 fs(points.at(i), 0.4, 0.2);
    // 	    s.search(std::inserter(neighbors, neighbors.end ()), fs);
    // 	    for (auto it : neighbors)
    // 	    {
    // 		indiceLocated.insert(boost::get<1>(it));
    // 	    }
    // 	}
    	
    // 	std::vector<int> pointLocated(indiceLocated.begin(), indiceLocated.end());
    // 	return pointLocated;
    // }
}
