#ifndef ALPHASHAPESURFACE_HPP
#define ALPAHSHAPESURFACE_HPP

#include "surface.hpp"

class AlphaShapeSurface : public Surface
{
public:
    void make(std::vector<Point>& points) override;
    float area() override;
    float volume() override;
    int locate(const std::vector<Point>& positions,
	       std::vector<int>& locatedPositions) override;
    
    void setAlphaValue(float alpha);

    
private:
    
    float alpha_ = 1.0;
    std::unique_ptr<Alpha_shape_3> dt_;
};

#endif
