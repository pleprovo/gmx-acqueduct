
#ifndef BRUTEFINDPAIRS_HPP
#define BRUTEFINDPAIRS_HPP

#include "find-pairs.hpp"

class BruteFindPairs : public FindPairs
{
public:
    std::vector<std::pair<int, int>> find(std::vector<Point>& points) override;
};

#endif
