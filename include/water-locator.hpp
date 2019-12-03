#ifndef WATERLOCATOR_HPP
#define WATERLOCATOR_HPP

class Point;

class WaterLocator
{
public:
    virtual void initialization(const std::vector<Point>& positions) = 0;
    virtual int locate(const std::vector<Point>& waterPositions,
		       std::vector<int>& results) = 0;
    
};

#endif
