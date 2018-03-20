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
				    
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
			      boost::no_property,
			      boost::property<boost::edge_index_t, std::size_t> > Graph;

// Descriptors
using vertex_descriptor = boost::graph_traits<Graph>::vertex_descriptor;
using edge_descriptor = boost::graph_traits<Graph>::edge_descriptor;

// Iterators
using vertex_iterator =  boost::graph_traits<Graph>::vertex_iterator;
using edge_iterator = boost::graph_traits<Graph>::edge_iterator;

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




#endif
