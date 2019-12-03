
#ifndef BRUTEFINDHYDROGENBONDS_HPP
#define BRUTEFINDHYDROGENBONDS_HPP

#include "find-hydrogen-bonds.hpp"

class BruteFindHydrogenBonds : public FindHydrogenBonds 
{
public:
    std::vector<HydrogenBond> find(const std::vector<Site>& sites) override;
    
};

#endif
