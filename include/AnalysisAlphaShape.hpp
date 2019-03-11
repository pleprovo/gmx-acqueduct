#ifndef STRATEGYALPHA_HPP
#define STRATEGYALPHA_HPP

#include "AnalysisInterface.hpp"

class StrategyAlpha : public AnalysisInterface
{
public:
    void initialize(Config config);
    Results execute(const Frame &frame);
};

#endif
