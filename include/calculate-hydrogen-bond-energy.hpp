#ifndef CALCULATEHYDROGENBONDENERGY_HPP
#define CALCULATEHYDROGENBONDENERGY_HPP

#include "utils.hpp"

class CalculateHydrogenBondEnergy
{
public:
    void setRadiusShift(float r_on, float r_off);
    void setAngleShift(float a_on, float a_off);
    // std::vector<HydrogenBond> calculate(const std::vector<std::pair<int, int>>& pairs,
    // 					const std::vector<Point>& sitePoints,
    // 					const std::vector<std::vector<Point_ptr>>& hydrogen,
    // 					const std::vector<SiteInfo_ptr>& sites);
    
private:
    float switch_radius(const float r);
    float switch_angle(const float a);
    
    float r_on_;
    float r_off_;
    float a_on_;
    float a_off_;
};

#endif
