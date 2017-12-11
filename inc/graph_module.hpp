/*
 *  Class for finding point inside the Alpha Shape Surface
 *  Pierre Leprovost, University of Oulu, Biocenter Oulu
 *  1/11/2015
 */

/*! \file
 * \brief
 * Declares FindElementInsideAlphaShape.
 *
 * \author Pierre Leprovost <leprovost.pierre@gmail.com>
 */

#ifndef GRAPHMODULE_HPP
#define GRAPHMODULE_HPP

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/range/iterator_range.hpp>


using Graph = boost::adjacency_list<boost::vecS,
				    boost::vecS,
				    boost::undirectedS,
				    boost::property<boost::vertex_index_t, int>,
				    boost::property<boost::edge_index_t, int> >;

// typedef subgraph < Graph > SubGraph;
using vertex_descriptor = boost::graph_traits<Graph>::vertex_descriptor;
using edge_descriptor = boost::graph_traits<Graph>::edge_descriptor;
using GraphTraits = boost::graph_traits<Graph>;

// Iterators
using vertex_iterator =  boost::graph_traits<Graph>::vertex_iterator;
using edge_iterator = boost::graph_traits<Graph>::edge_iterator;
using VertexIndexMap =  boost::property_map<Graph, boost::vertex_index_t>::type;
using EdgeIndexMap =  boost::property_map<Graph, boost::edge_index_t>::type;

using vertex_component_map = boost::shared_ptr<std::vector<unsigned long>>;

using ComponentGraph = boost::filtered_graph<Graph,
					     std::function<bool(Graph::edge_descriptor)>,
					     std::function<bool(Graph::vertex_descriptor)> >;



using Edge = std::pair<int, int>;

class GraphModule
{
public:
    GraphModule();
    void addGraph(const std::vector<Edge>& edgeList);
    
    
private:
    std::vector<Graph> graphFrames_;
    
};

#endif

std::vector<ComponentGraph> connected_components_subgraphs(const Graph &g);
float averageDegree(const Graph &g);
int writeGraphEdges(const Graph& g, std::stringstream& oss);




