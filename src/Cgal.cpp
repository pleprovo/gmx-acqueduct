
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
	
	std::set<int> indiceLocated;
	std::vector<cgal::Point_and_int> neighbors;
	for (unsigned  int i = 0; i < points.size(); i++)
	{
	    Fuzzy_sphere_3 fs(points.at(i), 0.4, 0.2);
	    s.search(std::inserter(neighbors, neighbors.end ()), fs);
	    for (auto it : neighbors)
	    {
		indiceLocated.insert(boost::get<1>(it));
	    }
	}
	// std::cout << neighbors.size() << std::endl;
	std::vector<int> pointLocated(indiceLocated.begin(), indiceLocated.end());
	return pointLocated;
    }

    
    /* Hydrogen Bond stuff */
    double switch_function(double r, double r_on, double r_off)
    {
	double sw = 0.0;

	if ( r_off > r && r > r_on ) {
	    sw = pow(pow(r,2)-pow(r_off,2),2);
	    sw *= pow(r_off,2)+2*pow(r,2)-3*pow(r_on,2);
	    sw /= pow(pow(r_off,2)-pow(r_on,2),3);	
	} else if ( r <= r_on ) {
	    sw = 1.0;
	}
	return sw;
    }


    double computeEnergy(const double r, const double cosine)
    {
	static float C = 3855; /* epsilon*sigma^6*sqrt(2/3) */
	static float D = 738; /* epsilon*sigma^4*sqrt(2/3) */
	static float theta_on = 0.25;
	static float theta_off = 0.0301;
	static float r_on = 5.5;
	static float r_off = 6.5;
	double r_switch = switch_function(r, r_on, r_off);
	double theta_switch = 1-switch_function(pow(cosine, 2), theta_off, theta_on); 
	double energy = -((C/pow(r, 6.0))-(D/pow(r, 4.0)));
	energy *= pow(cosine,4.0);
	energy *= r_switch;
	energy *= theta_switch;
	return energy;
    }

    
    double computeEnergy1(const double r, const double cosine,
			  const double r_on = 5.5,
			  const double r_off = 6.5,
			  const double theta_on = 0.25,
			  const double theta_off = 0.0301)
    {
	static float C = 3855; /* epsilon*sigma^6*sqrt(2/3) */
	static float D = 738; /* epsilon*sigma^4*sqrt(2/3) */
	// static float theta_on = 0.25;
	// static float theta_off = 0.0301;
	// static float r_on = 5.5;
	// static float r_off = 6.5;
	double r_switch = switch_function(r, r_on, r_off);
	double theta_switch = 1-switch_function(pow(cosine, 2), theta_off, theta_on); 
	double energy = -((C/pow(r, 6.0))-(D/pow(r, 4.0)));
	energy *= pow(cosine,4.0);
	energy *= r_switch;
	energy *= theta_switch;
	return energy;
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
		    energies.at(i) = computeEnergy1(distance, coss.at(i));
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
	    if ( e.length < 6.0 ) {
		edges.push_back(e);
	    }
	}

	
	return edges;
    }
}

