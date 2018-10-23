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
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/graph_traits.hpp>

#include <boost/range/iterator_range.hpp>

#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/successive_shortest_path_nonnegative_weights.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_utility.hpp>

#include <fstream>
#include <iostream>

using Traits = boost::adjacency_list_traits<boost::listS, boost::vecS, boost::directedS>;

struct Atom
{
    int id;
    Traits::vertex_descriptor predecessor;
    double distance;
    boost::default_color_type color;
};

struct HydrogenBond
{
    double length;
    double angle;
    double energy;
    double residual_energy;
    Traits::edge_descriptor reverse_edge;
};

using Graph = boost::adjacency_list<boost::listS, boost::vecS, boost::directedS,
				      Atom, HydrogenBond>;

class GraphModule
{
public:
    GraphModule();
    GraphModule(const int n);
    void add_vertex(const Atom &v);
    void add_edge(const int, const int, const HydrogenBond &e);
    double max_flow(const int, const int);
    void clear();
private:
    Graph g_;
};

#endif
