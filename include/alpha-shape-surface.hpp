#ifndef ALPHASHAPESURFACE_HPP
#define ALPAHSHAPESURFACE_HPP

#include "surface.hpp"

class AlphaShapeSurface : public Surface
{
public:
    virtual void make(std::vector<Point>& points) override;
    virtual float area() override;
    virtual float volume() override;
    virtual int locate(const std::vector<Point>& positions,
		       std::vector<int>& locatedPositions) override;
    
    void setAlphaValue(float alpha);

    
private:
    
    float alpha_ = 1.0;
    std::unique_ptr<Alpha_shape_3> dt_;
};

#endif
