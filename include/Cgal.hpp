
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
#include <functional>

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

    struct Site;
    
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
     * Switch function 
     */
    double switch_function(const double r, const double r_on, const double r_off);

    
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

    struct add_x {
	add_x(int x) : x(x) {}
	int operator()(int y) const { return x + y; }

    private:
	int x;
    };
    
     // Hydrogen Bond Energy Functor
    struct computeEnergy
    {
        computeEnergy(double r_on = 5.5,
		      double r_off = 6.5,
		      double theta_on = 0.25,
		      double theta_off = 0.0301) : r_on(r_on), 
						   r_off(r_off), 
						   theta_on(theta_on), 
						   theta_off(theta_off) {}
	
	double operator()(double r, double cosine) const
	    {
		static float C = 3855; /* epsilon*sigma^6*sqrt(2/3) */
		static float D = 738; /* epsilon*sigma^4*sqrt(2/3) */
		double r_switch = switch_function(r, r_on, r_off);
		double theta_switch = 1-switch_function(pow(cosine, 2),
							theta_off, theta_on); 
		double energy = -((C/pow(r, 6.0))-(D/pow(r, 4.0)));
		energy *= pow(cosine,4.0);
		energy *= r_switch;
		energy *= theta_switch;
		return energy;
	    }

    private:
	double r_on;
	double r_off;
	double theta_on ;
	double theta_off;
    };
    

    // Edge Energy Functor
    struct hydrogenBond
    {
	double length;
	double angle;
	double energy;
	
    };
    
    template<typename TEnergyFunctor, class Triangulation>
    struct computeEdgeEnergy
    {
	typedef typename Triangulation::edge_iterator Edge;
	
	computeEdgeEnergy(TEnergyFunctor &energyFunction)
	    : energyFunction(energyFunction) {}

	hydrogenBond operator() (const Edge &edge)
	    {
	        auto v1 = edge->first->vertex(edge->second);
		auto v2 = edge->first->vertex(edge->third);
		
		std::vector<double> energies(4, 0.0);
		auto OO = v1->point() - v2->point();
	    
		float distance = 10.0*CGAL::sqrt(OO*OO);

		auto OH11 = v1->info().h1 - v1->point();
		auto OH12 = v1->info().h2 - v1->point();
		auto OH21 = v2->info().h1 - v2->point();
		auto OH22 = v2->info().h2 - v2->point();
	    
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
			energies.at(i) = energyFunction(distance, coss.at(i));
		    }
		}
	    
		auto result = std::max_element(energies.begin(), energies.end());

		int ind = std::distance(energies.begin(), result);
		double energy = energies.at(ind);
		double angle = coss.at(ind);
		double length = distance;
		return hydrogenBond{length, angle, energy};

	    }
	
    private:
	TEnergyFunctor energyFunction;
	
    };
}
