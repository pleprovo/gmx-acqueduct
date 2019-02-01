
#ifndef WATERSEARCHINTERFACE_HPP
#define WATERSEARCHINTERFACE_HPP

class Point;

class WaterSearchInterface
{
public:
    std::vector<int> search(const std::vector<Point> &points) = 0;

}

#endif 
