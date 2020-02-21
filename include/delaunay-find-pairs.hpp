#ifndef DELAUNAYFINDPAIRS_HPP
#define DELAUNAYFINDPAIRS_HPP

#include "find-pairs.hpp"

class DelaunayFindPairs : public FindPairs
{
public:
    std::vector<std::pair<int, int>> find(std::vector<Point>& points) override;
};

#endif
