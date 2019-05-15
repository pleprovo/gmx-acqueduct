
#include <iostream>
#include <fstream>
#include <cmath>

#include "graph_module.hpp"
#include "Cgal.hpp"

/* Hydrogen Bond stuff */
double switch_function1(double r, double r_on, double r_off)
{
    double sw = 0.0;

    if ( r_off > r && r > r_on ) {
	sw = pow(pow(r,2)-pow(r_off,2),2);
	sw *= (pow(r_off,2)+2*pow(r,2)-3*pow(r_on,2));
	sw /= pow(pow(r_off,2)-pow(r_on,2),3);
	    
    } else if ( r <= r_on ) {
	sw = 1.0;
    }
    return sw;
}


double computeEnergy1(const double r, const double cosine)
{
    static float C = 3855; /* epsilon*sigma^6*sqrt(2/3) */
    static float D = 738; /* epsilon*sigma^4*sqrt(2/3) */
    static float theta_on = 0.25;
    static float theta_off = 0.0301;
    static float r_on = 5.5;
    static float r_off = 6.5;
    double r_switch = switch_function1(r, r_on, r_off);
    double theta_switch = 1-switch_function1(pow(cosine, 2), theta_off, theta_on); 
    double energy = -((C/pow(r, 6.0))-(D/pow(r, 4.0)));
    energy *= pow(cosine, 4.0);
    energy *= r_switch;
    energy *= theta_switch;	
    return energy;
}

int main (int argc, char *argv[])
{
    std::cout << " --- TESTS --- " << std::endl;
    std::cout << " --- Test Distance Switch : ";
    
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
    std::cout << " DONE --- " << std::endl;
    
    std::cout << " --- Test Angle Switch : ";
    
    // Testing angle switch
    std::ofstream oss2;
    oss2.open("switch_angle.dat");
    for (double cos = 0.0; cos < 1.0; cos+=0.001) {
    	oss2 << cos << " ";	
    }
    oss2 << "\n";
    for (double cos = 0.0; cos < 1.0; cos+=0.001) {
    	oss2 << 1-switch_function1(cos, 0.0301, 0.25) << " ";	
    }
    oss2.close();
    std::cout << " DONE --- " << std::endl;
    
    // Testing Hydrogen Bond potential
    std::cout << " --- Test Hydrogen Bond Potential : ";
    std::ofstream oss;
    oss.open("potential.dat");
    for (double r = 2.5; r < 7.0; r+=0.001) {
    	for (double theta = 0.0; theta < 1.0; theta+=0.01) {
    	    oss << computeEnergy1(r, theta) << " " ;
    	}
        oss << "\n";
    }
    oss.close();
    std::cout << " DONE --- " << std::endl;

    
    // Testing water water hb
    std::cout << " --- Test Hydrogen Bond Potential with CGAL : ";
    cgal::Point_3 O1(0.0, 0.0, 0.0);
    cgal::Point_3 H11(1.0, 0.0, 0.0);
    cgal::Point_3 H12(-0.5, 0.866, 0.0);
    cgal::Point_3 O2(3.5, 0.0, 0.0);
    cgal::Point_3 H21(4.5, 0.86, 0.0);
    cgal::Point_3 H22(4.5, -0.86, 0.0);

    using Vector_3 = CGAL::Vector_3<cgal::K>;
    
    Vector_3 O1O2 = O2 - O1;
    Vector_3 O2O1 = -O1O2;
    
    double distance = CGAL::sqrt(O1O2*O1O2);
    
    Vector_3 OH11 = H11 - O1;
    Vector_3 OH12 = H12 - O1;
    Vector_3 OH21 = H21 - O2;
    Vector_3 OH22 = H22 - O2;

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

    std::vector<std::pair<int, double> > energies;
    for (unsigned int i = 0; i< coss.size(); i++) {
    	if (coss.at(i) > 0.0) {
    	    energies.push_back(std::pair<int,double>(i,
						     computeEnergy1(distance,
								    coss.at(i))));
    	}
    }
    
    std::cout << " DONE --- " << std::endl;
    
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
