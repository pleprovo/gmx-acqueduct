
#include "calculate-hydrogen-bond-energy.hpp"

#include <math.h>

void CalculateHydrogenBondEnergy::setRadiusShift(float r_on, float r_off)
{
    r_on_ = r_on*r_on;
    r_off_ = r_off*r_off;
}

void CalculateHydrogenBondEnergy::setAngleShift(float a_on, float a_off)
{
    a_on_ = a_on*a_on;
    a_off_ = a_off*a_off;
}

// float CalculateHydrogenBondEnergy::calculate(const HydrogenBond& hb)
// {
//     static float C = 3855; /* epsilon*sigma^6*sqrt(2/3) */
//     static float D = 738; /* epsilon*sigma^4*sqrt(2/3) */
//     const float& r = hb.length;
//     const float& a = hb.angle;
//     double energy = ((C/pow(r, 6.0))-(D/pow(r, 4.0)));
//     energy *= pow(a, 4.0);
//     energy *= switch_radius(r*r);
//     energy *= switch_angle(a*a);
//     return energy;
// }
    

float CalculateHydrogenBondEnergy::switch_radius(const float r)
{
    double sw = 0.0;
    	
    if ( r_off_ > r && r > r_on_ ) {
	sw = pow(pow(r,2)-pow(r_off_,2),2);
	sw *= pow(r_off_,2)+2*pow(r,2)-3*pow(r_on_,2);
	sw /= pow(pow(r_off_,2)-pow(r_on_,2),3);	
    } else if ( r < r_on_ ) {
	sw = 1.0;
    }
    return sw;
}


float CalculateHydrogenBondEnergy::switch_angle(const float a)
{
    double sw = 0.0;
		
    if ( a_off_ < a && a < a_on_ ) {
	sw = pow(pow(a,2)-pow(a_off_,2),2);
	sw *= pow(a_off_,2)+2*pow(a,2)-3*pow(a_on_,2);
	sw /= pow(pow(a_off_,2)-pow(a_on_,2),3);
	sw *= -1.0;
    } else if ( a > a_on_ ) {
	sw = 1.0;
    }
    return sw;
}
