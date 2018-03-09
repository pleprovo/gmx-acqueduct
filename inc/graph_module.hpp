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

#include <boost/range/iterator_range.hpp>
/*
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/graph_utility.hpp> 
*/

#include <fstream>
#include <iostream>

using Traits = boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::undirectedS>;

				    
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
			      boost::no_property,
			      boost::property<boost::edge_index_t, std::size_t> > Graph;

// Descriptors
using vertex_descriptor = boost::graph_traits<Graph>::vertex_descriptor;
using edge_descriptor = boost::graph_traits<Graph>::edge_descriptor;

// Iterators
using vertex_iterator =  boost::graph_traits<Graph>::vertex_iterator;
using edge_iterator = boost::graph_traits<Graph>::edge_iterator;

using HydrogenBond = std::pair<std::pair<int, int>, float>;

/*
typedef subgraph < Graph > SubGraph;
using vertex_component_map = boost::shared_ptr<std::vector<unsigned long>>;

using ComponentGraph = boost::filtered_graph<Graph,
					     std::function<bool(Graph::edge_descriptor)>,
					     std::function<bool(Graph::vertex_descriptor)> >;

std::vector<ComponentGraph> connected_components_subgraphs(const Graph &g);

Subgraph buildSubgraph(const std::vector<int> &vertices, const Subgraph &g);
*/
/*
template <typename T>
float averageDegree(T &g) {
    vertex_iterator vi, vi_end, next;
    boost::tie(vi, vi_end) = boost::vertices(g);
    float average = 0.0;
    for (next = vi; vi != vi_end; vi = next) {
	++next;
        average += boost::degree(*vi, g);
    }
    return average/boost::num_vertices(g);

}
*/

void AddBidirectionalEdge(Graph& graph,
			  unsigned int source,
			  unsigned int target,
			  float weight,
                          std::vector<edge_descriptor>& reverseEdges,
			  std::vector<float>& capacity)
{
    // Add edges between grid vertices. We have to create the edge and the reverse edge,
    // then add the reverseEdge as the corresponding reverse edge to 'edge', and then add 'edge'
    // as the corresponding reverse edge to 'reverseEdge'
    int nextEdgeId = num_edges(graph);

    edge_descriptor edge;
    bool inserted;

    boost::tie(edge,inserted) = add_edge(source, target, nextEdgeId, graph);
    if(!inserted)
    {
        std::cerr << "Not inserted!" << std::endl;
    }
    edge_descriptor reverseEdge = add_edge(target, source, nextEdgeId + 1, graph).first;
    reverseEdges.push_back(reverseEdge);
    reverseEdges.push_back(edge);
    capacity.push_back(weight);

    // Not sure what to do about reverse edge weights
    capacity.push_back(weight);
//    capacity.push_back(0);
}

#endif
