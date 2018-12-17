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

#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_utility.hpp>

#include <boost/graph/connected_components.hpp>
#include <boost/graph/subgraph.hpp>

#include <fstream>
#include <iostream>

using Traits = boost::adjacency_list_traits<boost::listS, boost::vecS, boost::directedS>;
using Utraits = boost::adjacency_list_traits<boost::listS, boost::vecS, boost::undirectedS>;

struct Atom
{
    int id;
    int atom_id;
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

using Ugraph = boost::adjacency_list<boost::listS, boost::vecS, boost::directedS,
				     Atom, HydrogenBond>;

using SubGraph = boost::subgraph<Graph>;

class MyVisitor : public boost::default_dfs_visitor
{
  public:
    MyVisitor() : vv(new std::vector<int>()) {}

    void discover_vertex(int v, const Graph &g) const {
        vv->push_back(g[v].id);
	std::cout << " Discovered : " << g[v].id << std::endl;
    }

    std::vector<int> &GetVector() const { return *vv; }

  private:
    std::vector<int> *vv;
};

class GraphModule
{
public:
    GraphModule();
    void add_vertex(const Atom &v);
    void add_edge(const int, const int, const HydrogenBond &e, const std::string mode = "full");
    void add_edge1(const int, const int, const HydrogenBond &e);
    double max_flow(const int, const int);
    std::vector<int> dfs(const int);
    void clear();
private:
    Graph g_;
    MyVisitor vis_;
};

#endif
