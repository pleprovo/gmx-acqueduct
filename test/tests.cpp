
#include <iostream>

#include "graph_module.hpp"
#include "alpha_shape_module.hpp" 

/* Hydrogen Bond stuff */
double switch_function1(double r, double r_on, double r_off)
{
    double sw = 0.0;

    if ( r_off >= r && r >= r_on ) {
	sw = (pow(pow(r,2)-pow(r_off,2),2)*(pow(r_off,2)+2*pow(r,2)-3*pow(r_on,2)))/pow(pow(r_off,2)-pow(r_on,2),3);
    } else if ( r < r_on ) {
	sw = 1.0;
    }
    return sw;
}


double computeEnergy1(const double r, const double cosine)
{
    static float C = 3855; /* epsilon*sigma^6*sqrt(2/3) */
    static float D = 738; /* epsilon*sigma^4*sqrt(2/3) */
    static float theta_on = 0.25;
    static float theta_off = 0.003;
    static float r_on = 5.5;
    static float r_off = 6.5;
    double energy = 0.0;
    double r_switch = switch_function1(r, r_on, r_off);
    double theta_switch = 1-switch_function1(pow(cosine, 2), theta_off, theta_on); 
    energy = -((C/pow(r, 6.0))-(D/pow(r, 4.0)))*r_switch*theta_switch*pow(cosine, 4.0);
    return energy;
}

int main (int argc, char *argv[])
{
    std::cout << " --- TESTS --- " << std::endl;

    Graph g;
    GraphModule gm;
    
    gm.add_vertex(Atom{0});
    gm.add_vertex(Atom{1});
    gm.add_vertex(Atom{2});
    gm.add_vertex(Atom{3});
    gm.add_vertex(Atom{4});
    gm.add_vertex(Atom{5});

    gm.add_edge(0, 1, HydrogenBond{3.5, 0.0, 4.0});
    gm.add_edge(1, 2, HydrogenBond{3.5, 0.0, 4.0});
    gm.add_edge(1, 3, HydrogenBond{3.5, 0.0, 4.0});
    gm.add_edge(3, 4, HydrogenBond{3.5, 0.0, 4.0});
    gm.add_edge(2, 4, HydrogenBond{3.5, 0.0, 4.0});
    gm.add_edge(4, 5, HydrogenBond{3.5, 0.0, 4.0});

    // std::cout <<  " Test Graph : " << std::endl;
    // boost::print_graph(g, boost::get(&Atom::resid, g));
    // std::cout << std::endl;

    // boost::dijkstra_shortest_paths(g, v0, boost::weight_map(get(&HydrogenBond::length,g)).distance_map(get(&Atom::distance, g)).predecessor_map(get(&Atom::predecessor, g)));
    
    // print_predecessor_path(g, v5);
    

    double flow1 = gm.max_flow(0, 5);
    std::cout << "Flow 1: " << flow1 << "\n";
    double flow2 = gm.max_flow(5, 0);
    std::cout << "Flow 2: " << flow2 << "\n";
    gm.clear();

    // Testing radius switch
    std::ofstream oss1;
    oss1.open("switch_radius.dat");
    for (double r = 2.0; r < 7.0; r+=0.1) {
	oss1 << r << " ";	
    }
    oss1 << "\n";
    for (double r = 2.0; r < 7.0; r+=0.1) {
	oss1 << switch_function1(r, 5.5, 6.5) << " ";
    }
    oss1.close();
	
    // Testing angle switch
    std::ofstream oss2;
    oss2.open("switch_angle.dat");
    for (double cos = 0.0; cos < 1.0; cos+=0.001) {
	oss2 << cos << " ";	
    }
    oss2 << "\n";
    for (double cos = 0.0; cos < 1.0; cos+=0.001) {
	oss2 << 1-switch_function1(cos, 0.25, 0.003) << " ";	
    }
    oss2.close();
    
    // Testing Hydrogen Bond potential
    std::ofstream oss;
    oss.open("potential.dat");
    for (double r = 2.0; r < 6.5; r+=0.1) {
	for (double theta = 0.0; theta < 1.0; theta+=0.1) {
	    oss << computeEnergy1(r, theta) << " " ;
	}
        oss << "\n";
    }
    oss.close();

    // Testing water water hb
    Point_3 O1(0.0, 0.0, 0.0);
    Point_3 H11(1.0, 0.0, 0.0);
    Point_3 H12(0.0, 0.1, 0.0);
    Point_3 O2(3.0, 0.0, 0.0);
    Point_3 H21(2.5, 0.86, 0.0);
    Point_3 H22(2.5, 0.86, 0.0);
    
    CGAL::Vector_3<K> O1O2 = O2 - O1;
    CGAL::Vector_3<K> O2O1 = O1 - O2;
    std::cout << O1O2 << " / " << O2O1 << std::endl;
    
    double distance = CGAL::sqrt(O1O2*O1O2);
    std::cout << "Distance : " << distance << std::endl;
    
    CGAL::Vector_3<K> OH11 = H11 - O1;
    CGAL::Vector_3<K> OH12 = H12 - O1;
    CGAL::Vector_3<K> OH21 = H21 - O2;
    CGAL::Vector_3<K> OH22 = H22 - O2;

    O1O2 = O1O2 / CGAL::sqrt(O1O2*O1O2);
    O2O1 = O2O1 / CGAL::sqrt(O2O1*O2O1);
    
    OH11 = OH11 / CGAL::sqrt(OH11*OH11);
    OH12 = OH12 / CGAL::sqrt(OH12*OH12);
    OH21 = OH21 / CGAL::sqrt(OH21*OH21);
    OH22 = OH22 / CGAL::sqrt(OH22*OH22);

    std::vector<double> coss;
    coss.push_back(OH11 * O1O2);
    coss.push_back(OH12 * O1O2);
    coss.push_back(OH21 * O2O1);
    coss.push_back(OH22 * O2O1);
    
    std::cout << "cos11 : " << coss.at(0) << ", "
	      << "cos12 : " << coss.at(1) << ", "
	      << "cos21 : " << coss.at(2) << ", "
	      << "cos22 : " << coss.at(3) << std::endl;

    std::vector<std::pair<int, double> > energies;
    for (unsigned int i = 0; i< coss.size(); i++) {
	if (coss.at(i) > 0.0) {
	    energies.push_back(std::pair<int,double>(i, computeEnergy1(distance, coss.at(i))));
	}
    }
    if (!energies.empty()) {
	for (auto& ener : energies) {
	    std::cout << ener.first << " : " << ener.second << std::endl; 
	}
	std::cout << std::endl;
    } else {
	std::cout << " Not HB at all ! " << std::endl;
    }
    // Testing FFTW
    // int N = data.size();
    // double in[N];
    // fftw_complex *out;
    // fftw_plan fft, ifft;

    // out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
    // fft = fftw_plan_dft_r2c_1d(N, in, out, FFTW_FORWARD);
    
    // fftw_execute(my_plan);    
    // fftw_destroy_plan(my_plan);
    
    // fftw_free(in);    
    // fftw_free(out);
    
    return 0;
}
