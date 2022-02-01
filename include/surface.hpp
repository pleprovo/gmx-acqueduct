#ifndef SURFACE_HPP
#define SURFACE_HPP

#include "cgal.hpp"


class Surface
{
public:
    virtual void make(std::vector<Point>& points) = 0;
    virtual float area() = 0;
    virtual float volume() = 0;
    virtual int locate(const std::vector<Point>& positions,
		       std::vector<int>& locatedPositions) = 0;
};

#endif
