#ifndef DELAUNAYFINDHYDROGENBONDS_HPP
#define DELAUNAYFINDHYDROGENBONDS_HPP

#include "find-hydrogen-bonds.hpp"

class DelaunayFindHydrogenBonds : public FindHydrogenBonds 
{
public:
    std::vector<HydrogenBond> find(const std::vector<Site>& sites) override;

};

#endif
