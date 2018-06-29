/*
 *  Class for finding point inside the Alpha Shape Surface
 *  Pierre Leprovost, University of Oulu, Biocenter Oulu
 *  1/11/2015
 */

/*! \file
 * \brief
 * Declares FindElementInsideAlphaShape.
 *
 * \author Pierre Leprovost <leprovost.pierre@gmail.com>
 */

#ifndef ALPHASHAPEMODULE_HPP
#define ALPHASHAPEMODULE_HPP
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>


// CGAL Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// typedef for Alpha Shape
typedef CGAL::Alpha_shape_vertex_base_3<K>                  Vb;
typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vb>     Vbh;
typedef CGAL::Alpha_shape_cell_base_3<K>                    Fb;
typedef CGAL::Triangulation_data_structure_3<Vbh,Fb>        Tds;
typedef CGAL::Delaunay_triangulation_3<K,Tds>               Delaunay;
typedef CGAL::Triangulation_hierarchy_3<Delaunay>           Delaunay_hierarchy;
typedef CGAL::Alpha_shape_3<Delaunay_hierarchy>             Alpha_shape_3;
typedef Alpha_shape_3::Facet		                    Facet;
typedef Alpha_shape_3::Classification_type                  Classification;
typedef K::Point_3					    Point_3;

/*!
  * Class that implement the method to find the water inside the protein.
  * The method depend on the library CGAL to calculated the surface, 
  * the surface in define as an alpha Shape. The elements inside the alpha shape
  * are also fiond with CGAL methods.
  */


class AlphaShapeModule
{
private:
    std::shared_ptr<Alpha_shape_3> alphaShape_;
        
public:   
    AlphaShapeModule();
    AlphaShapeModule(const std::vector<Point_3> &pointVector, const float &alpha);
    virtual ~AlphaShapeModule();

    // Modifier
    void build(const std::vector<Point_3> &pointVector, const float &alpha);
    void update(const std::vector<Point_3> &pointVector);

    // Special functions
    std::vector<int> locate(const std::vector<Point_3> &pointVector, const bool &location) const;
  
    // Accessors 
    float volume();    
    void writeOff(std::string &outputString);
    
};

#endif
