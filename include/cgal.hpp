
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
#include <CGAL/boost/iterator/counting_iterator.hpp>

#include <tuple>

// #ifdef CGAL_LINKED_WITH_TBB
// typedef CGAL::Parallel_tag Concurrency_tag;
// #else
// typedef CGAL::Sequential_tag Concurrency_tag;
// #endif

// #include "utils.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector_3;

using Point_ptr = std::shared_ptr<Point>;

// Definition for Alpha Shape 3
typedef CGAL::Triangulation_data_structure_3<
    CGAL::Triangulation_hierarchy_vertex_base_3<CGAL::Alpha_shape_vertex_base_3<K>>,
    CGAL::Alpha_shape_cell_base_3<K>,
    CGAL::Parallel_tag> Tdsa;
typedef CGAL::Delaunay_triangulation_3<K,
				       Tdsa,
				       CGAL::Fast_location> Delaunay;
typedef CGAL::Alpha_shape_3<Delaunay> Alpha_shape_3;

// Definition for Delaunay with info
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Parallel_tag> Tds;
typedef CGAL::Delaunay_triangulation_3<K,
				       Tds,
				       CGAL::Parallel_tag> DelaunayWithInfo;

// Adapt spatial sort for delaunay with info
typedef CGAL::Spatial_sort_traits_adapter_3<K,
					    CGAL::Pointer_property_map<Point>::type > Search_traits_3;


