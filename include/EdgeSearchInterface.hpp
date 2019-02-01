

#ifndef EDGESEARCHINTERFACE_HPP
#define EDGESEARCHINTERFACE_HPP

class Edge;
class Node;

class EdgeSearchInterface
{
public:
    std::vector<Edge> search(const std::vector<Node> nodes) = 0;
    
};


#endif
