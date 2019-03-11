#ifndef ANALYSISINTERFACE_HPP
#define ANALYSISINTERFACE_HPP

#include <gromacs/trajectoryanalysis.h>

struct Config {
    int plop;
};

struct Frame {
    const gmx::Selection &waters;
    const gmx::Selection &points;
    const gmx::Selection &proteins; 
};

struct Results {
    int numVertices;
    int numEdges;
};

template <class T>
std::vector<T> fromGmxtoCgalPosition(const gmx::ConstArrayRef<rvec> &coordinates,
				     const int increment=1)
{
    std::vector<T> cgalPositionVector;   
    for (unsigned int i = 0; i < coordinates.size(); i += increment)
    {
	cgalPositionVector.push_back(T(coordinates.at(i)[XX],
				       coordinates.at(i)[YY],
				       coordinates.at(i)[ZZ]));
    }
    return cgalPositionVector;
}

class AnalysisInterface
{
public:
    virtual void initialize(Config config) = 0;
    virtual Results execute(const Frame &frame) = 0;
};

#endif
/*

Alpha Shape needs double alpha, vector alpha, vector points
KDTree needs double cutoff, vector alpha, vector search

K Neighbor needs int K, vector points
Delaunay needs vector points 

 */
