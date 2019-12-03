#ifndef ALPHASHAPEWATERLOCATOR_HPP
#define ALPAHSHAPEWATERLOCATOR_HPP

#include "water-locator.hpp"

class AlphaShapeWaterLocator : public AlphaShapeWaterLocator
{
public:
    void initialization(const std::vector<Point>& positions) override;
    int locate(const std::vector<Point>& waterPositions,
	       std::vector<int>& results) override;

    
};

#endif
