#ifndef CALCULATEHYDROGENBONDENERGY_HPP
#define CALCULATEHYDROGENBONDENERGY_HPP

#include "find-hydrogen-bonds.hpp"

class CaculateHydrogenBondEnergy
{
public:
    float calculate(const HydrogenBond& hb);
private:
    
    distanceCuton_;
    distanceCutoff_;
    angleCuton_;
    angleCutoff_;
}
