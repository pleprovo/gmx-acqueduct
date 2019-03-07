
#include "Cgal.hpp"

namespace cgal {

    float getVolume(Alpha_shape_3 &s)
    {
	float volume = 0;
	
	for(Alpha_shape_3::Finite_cells_iterator it = s.finite_cells_begin(); 
	    it != s.finite_cells_end(); it++)
	{
	    if(s.classify(it)==Alpha_shape_3::INTERIOR)
	    { 
		volume += s.tetrahedron(it).volume();
	    }
	}
	
	return volume;
    }
    
    std::vector<int> searchPoints(Alpha_shape_3 &s, std::vector<Point_3> points)
    {
	std::vector<int> pointLocated;      
    
	for (unsigned int i = 0; i < points.size(); i++)
	{
	    if (s.classify(points.at(i)) == Alpha_shape_3::INTERIOR)
	    {
		pointLocated.push_back(i);
	    }
	}

	return pointLocated;
    }

    std::vector<int> searchPoints(Kd_tree_3 &s, std::vector<Point_3> points)
    {
	std::vector<int> pointLocated;

	for (unsigned  int i = 0; i < points.size(); i++)
	{
	    Fuzzy_sphere_3 fs(points.at(i), 3.0, 1.0);
	    //s.search()
	}

	return pointLocated;
    }
    
    std::vector<edge> analyseEdges(DelaunayWithInfo &g)
    {
	std::vector<edge> edges;
	CGAL::Vector_3<K> OO, OH11, OH12, OH21, OH22, OH1, OH2;
	DelaunayWithInfo::Vertex_handle v1, v2;
	for(DelaunayWithInfo::Finite_edges_iterator ei=g.finite_edges_begin();
	    ei!=g.finite_edges_end(); ++ei) {
	    v1 = ei->first->vertex(ei->second);
	    v2 = ei->first->vertex(ei->third);

	    edge e;

	    std::vector<double> energies(4, 0.0);
	    OO = v1->point() - v2->point();
	    
	    float distance = 10.0*CGAL::sqrt(OO*OO);

	    OH11 = v1->info().h1 - v1->point();
	    OH12 = v1->info().h2 - v1->point();
	    OH21 = v2->info().h1 - v2->point();
	    OH22 = v2->info().h2 - v2->point();
	    
	    OO = OO / CGAL::sqrt(OO*OO);
	    OH11 = OH11 / CGAL::sqrt(OH11*OH11);
	    OH12 = OH12 / CGAL::sqrt(OH12*OH12);
	    OH21 = OH21 / CGAL::sqrt(OH21*OH21);
	    OH22 = OH22 / CGAL::sqrt(OH22*OH22);

	    std::vector<double> coss;
	    coss.push_back(OH11 * OO);
	    coss.push_back(OH12 * OO);
	    coss.push_back(OH21 * -OO);
	    coss.push_back(OH22 * -OO);
	    
	    for (unsigned int i = 0; i < coss.size(); i++) {
		if (coss.at(i) > 0.0) {
		    // energies.at(i) = computeEnergy(distance, coss.at(i));
		    // std::invoke
		}
	    }
	    
	    auto result = std::max_element(energies.begin(), energies.end());

	    int ind = std::distance(energies.begin(), result);
	    e.energy = energies.at(ind);
	    e.energy = 0.0;
	    e.angle = coss.at(ind);
	    e.length = distance;
	    
	    if ( ind < 2 ) {
		e.i = v1->info().id;
		e.j = v2->info().id;
	    }
	    else
	    {
		e.i = v2->info().id;
		e.j = v1->info().id;
	    }
	    if ( e.length < 3.5 ) {
		edges.push_back(e);
	    }
	}

	
	return edges;
    }
}

