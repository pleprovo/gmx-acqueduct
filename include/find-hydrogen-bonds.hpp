#ifndef FINDHYDROGENBONDS_HPP
#define FINDHYDROGENBONDS_HPP

#include <tuple>

#include "site.hpp"

struct HydrogenBond {
    std::pair<Site_ptr, Site_ptr> sites;
    float length;
    float angle;
    float energy;
};

class FindHydrogenBonds
{
public:
    void setDistanceCutoff(const float distance = 6.5)
	{
	    distanceCutoff_ = distance;
	}

    void setAngleCutoff(const float angle = 30 )
	{
	    angleCutoff_ = angle;
	}
    
    virtual std::vector<HydrogenBond> find(const std::vector<Site>& sites) = 0;
    
protected:
    float distanceCutoff_;
    float angleCutoff_;
    
};

#endif
