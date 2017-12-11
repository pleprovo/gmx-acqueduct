#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/periodic_3_triangulation_3_io.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K>       Gt;
typedef CGAL::Periodic_3_Delaunay_triangulation_3<Gt>             P3DT3;
typedef P3DT3::Point             Point;
typedef P3DT3::Iso_cuboid        Iso_cuboid;
typedef P3DT3::Vertex_handle     Vertex_handle;
typedef P3DT3::Cell_handle       Cell_handle;
typedef P3DT3::Locate_type       Locate_type;
typedef K::Point_3               Point_3;

#ifndef PERIODICTRIANGULATION_HPP
#define PERIODICTRIANGULATION_HPP


class PeriodicTriangulation : public P3DT3
{
public:
    PeriodicTriangulation(std::vector<Point>&, Iso_cuboid&);
    void set_domain(float, float, float, float, float, float);
    Iso_cuboid domain();
    std::vector<std::pair<int, int> > filter_edges(float);
    
private:
    Iso_cuboid domain_;
    
};

#endif
