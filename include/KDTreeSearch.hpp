
#ifndef KDTREESEARCH_HPP
#define KDTREESEARCH_HPP

#include "WaterSearchInterface.hpp"

class KDTreeSearch : public WaterSearchInterface
{
private:
    std::shared_ptr<KDTree> kdTree_;
    std::vector<Point> 
public:
    
    std::vector<int> search(const std::vector<Point> &points);

}


#endif 
