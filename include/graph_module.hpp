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
    int resid;
    std::string name;
    bool isDonor;
    bool isWater;
    Traits::vertex_descriptor predecessor;
    double distance;
    boost::default_color_type color;
};

struct HydrogenBond
{
    int id;
    double length;
    double angle;
    double energy;
    double residual_energy;
    Traits::edge_descriptor reverse_edge;
};

using Graph = boost::adjacency_list<boost::listS, boost::vecS, boost::directedS,
				      Atom, HydrogenBond>;

template <typename S>
void add_bidirectional_edge(int u, int v, S s, Graph &g)
{
    auto e1 = boost::add_edge(u, v, s, g).first;
    auto e2 = boost::add_edge(v, u, s, g).first;
    g[e1].reverse_edge = e2;
    g[e2].reverse_edge = e1;
}

void print_predecessor_path(Graph &g, Traits::vertex_descriptor v);
double do_max_flow(Graph &g, Graph::vertex_descriptor &source, Graph::vertex_descriptor &sink);

    
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
