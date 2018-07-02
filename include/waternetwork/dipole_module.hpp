



#ifndef DIPOLEMODULE_HPP
#define DIPOLEMODULE_HPP

#include "gromacs/math/vec.h"
#include "gromacs/utility/arrayref.h"
#include <vector>

class DipoleModule
{
public:
    DipoleModule();
    void initialise(const gmx::ConstArrayRef<rvec>& waterCoordinates);
    void analyseFrame(const gmx::ConstArrayRef<rvec>& waterCoordinates);
    real averageFrame(const std::vector<int>&);
    std::vector<real> average(const std::vector<int>&);
    
private:
    std::vector<rvec> initDipoles_;
    std::vector<std::vector<real> > dipolesFrames_;
};

#endif 
