#ifndef FINDPAIRS_HPP
#define FINDPAIRS_HPP

#include "cgal.hpp"

#include <vector>

class FindPairs
{
public:  
    virtual std::vector<std::pair<int, int>> find(std::vector<Point>& points) = 0;

};

#endif
