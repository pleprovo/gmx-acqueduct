#ifndef ANALYSISNEIGHBORS_HPP
#define ANALYSISNEIGHBORS_HPP

#include "AnalysisInterface.hpp"

class AnalysisNeighbors : public AnalysisInterface
{
public:
    void initialize(Config config);
    void execute(const Frame &frame);
};

#endif
