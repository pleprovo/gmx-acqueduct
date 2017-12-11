
#include "periodictriangulation.hpp"

PeriodicTriangulation::PeriodicTriangulation(std::vector<Point>& points, Iso_cuboid& box)
    : P3DT3(points.begin(), points.end(), box) {

}

void PeriodicTriangulation::set_domain(float x, float y, float z, float u, float v, float w) {
    this->domain_ = Iso_cuboid(x, y, z, u, v, w);
}

Iso_cuboid PeriodicTriangulation::domain() {
    return this->domain_;
}

std::vector<std::pair<int, int> > PeriodicTriangulation::filter_edges(float cutOff) {
    std::vector<std::pair<int, int> > edges;

    for (P3DT3::Edge_iterator ei = this->edges_begin(); ei != this->edges_end(); ei++) { 
	int i1 = (ei->second + 1) % 3;
	int i2 = (ei->second + 2) % 3;
	P3DT3::Vertex_handle v1 = ei->first->vertex(i1);
	P3DT3::Vertex_handle v2 = ei->first->vertex(i2);

	float l = CGAL::squared_distance(v1->point(), v2->point());
	
	if ( l < cutOff ) {
	    edges.push_back(std::pair<int, int> (i1, i2));
	}
	
    }
    
    return edges;
}

