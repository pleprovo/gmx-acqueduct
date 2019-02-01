#ifndef ALPHASHAPESEARCH_HPP
#define ALPHASHAPESEARCH_HPP

#include "Cgal.hpp"
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>


typedef CGAL::Alpha_shape_vertex_base_3<K>                  Vb;
typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vb>     Vbh;
typedef CGAL::Alpha_shape_cell_base_3<K>                    Fb;
typedef CGAL::Triangulation_data_structure_3<Vbh,Fb>        Tds;
typedef CGAL::Delaunay_triangulation_3<K,Tds, CGAL::Fast_location> Delaunay;
typedef CGAL::Alpha_shape_3<Delaunay>                       Alpha_shape_3;


class AlphaShapeSearch : public WaterSearchInterface
{
private:
    /*
     *  CGAL Alpha Shape Obkject Pointer  
     */
    std::shared_ptr<Alpha_shape_3> alphaShape_;

    /*
     * Alpha value for alpha shape computation
     */
    float alpha_;

    /*
     * Location selector for points search 
     */
    bool locator_;
        
public:   
    AlphaShapeSearch();
    virtual ~AlphaShapeSearch();
    /*
     *  Set the alpha value for the alpha shape computation
     */
    void setAlpha(const float alpha);

    /*
     * Set locator interior of exterior search 
     */
    void setLocator(const bool locator);

    /*
     * Compute and return the volume of the alpha shape (zero is no alpha shape)
     */
    float getVolume();

    /*
     * Compute the alpha shape with the given set of point and alpha value
     */
    void build(const std::vector<Point_3> &pointVector);

    /*
     * Search if any of the given points is located inside or outside the alpha shape
     */
    std::vector<int> search(const std::vector<Point> &points);

    /*
     * Write the alpha shape in OFF format to string
     */
    void writeOff(std::string &outputString);
    
};

#endif
