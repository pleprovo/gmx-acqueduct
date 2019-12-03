
#include "cgal.hpp"

#include <future>

namespace cgal {
    double switch_function_radius(const double r, const double r_on, const double r_off)
    {
	double sw = 0.0;
		
	if ( r_off > r && r > r_on ) {
	    sw = pow(pow(r,2)-pow(r_off,2),2);
	    sw *= pow(r_off,2)+2*pow(r,2)-3*pow(r_on,2);
	    sw /= pow(pow(r_off,2)-pow(r_on,2),3);	
	} else if ( r < r_on ) {
	    sw = 1.0;
	}
	return sw;
    }

    double switch_function_angle(const double a, const double a_on, const double a_off)
    {
	double sw = 0.0;
		
	if ( a_off < a && a < a_on ) {
	    sw = pow(pow(a,2)-pow(a_off,2),2);
	    sw *= pow(a_off,2)+2*pow(a,2)-3*pow(a_on,2);
	    sw /= pow(pow(a_off,2)-pow(a_on,2),3);
	    sw *= -1.0;
	} else if ( a > a_on ) {
	    sw = 1.0;
	}
	return sw;
    }

    double computeEnergy(const double r,
			 const double cosine,
			 const double r_on,
			 const double r_off,
			 const double theta_on,
			 const double theta_off)
    {
	static float C = 3855; /* epsilon*sigma^6*sqrt(2/3) */
	static float D = 738; /* epsilon*sigma^4*sqrt(2/3) */
	double radius_switch = switch_function_radius(r*r, r_on*r_on, r_off*r_off);
	double angle_switch = switch_function_angle(cosine*cosine,
						    theta_off, theta_on); 
	double energy = ((C/pow(r, 6.0))-(D/pow(r, 4.0)));
	energy *= pow(cosine, 4.0);
	energy *= radius_switch;
	energy *= angle_switch;
	return energy;
    }

    float getVolume(Alpha_shape_3 &s)
    {
	float volume = 0.0;
	
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

    std::vector<int> filterPoints (const std::vector<Point_3>& points,
    				   Alpha_shape_3& s,
    				   const int start,
    				   const int stop)
    {
    	std::vector<int> pointLocated;
    	for ( int i = start; i < stop; i++)
    	{
    	    if (s.classify(points.at(i)) == Alpha_shape_3::INTERIOR)
    	    {
    		pointLocated.push_back(i);
    	    }
    	}
    	return pointLocated;
    }

    int filterPointsParallel(const std::vector<Point_3>& points,
			     Alpha_shape_3& s,
			     std::vector<int>& results)
    {
	/* Make chunks */
	int n = std::thread::hardware_concurrency();
	int chunk_size = points.size() / n;
	/* Make Futures */
	std::vector< std::future< std::vector< int> > > futures;
	for ( int k = 0; k < n-1; k++ )
	{
	    futures.push_back(std::async(std::launch::async,
					 filterPoints, std::ref(points), std::ref(s),
					 chunk_size*k, chunk_size*k+chunk_size));
	}
	
	results = filterPoints(points, s, chunk_size*(n-1), points.size());
	
	/* Gather results */
	for ( std::future< std::vector<int> >& fut : futures )
	{
	    fut.wait();
	}

	for (  std::future< std::vector<int> >& fut : futures )
	{
	    std::vector<int> chunk = fut.get();
	    results.insert(results.end(),
			   std::make_move_iterator(chunk.begin()),
			   std::make_move_iterator(chunk.end()));
	}
	
	return results.size();

    }
}
