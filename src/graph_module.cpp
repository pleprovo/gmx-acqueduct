
#include "graph_module.hpp"

GraphModule::GraphModule()
{
    g_ = Graph();
}

void GraphModule::add_vertex(const Atom &v)
{
    boost::add_vertex(v, g_);
}


void GraphModule::add_edge(const int u, const int v, const HydrogenBond &e, const std::string mode)
{
    auto e1 = boost::add_edge(u, v, e, g_).first;
    auto er = e;
    if (mode == "null") {
        er.energy = 0.0;
    }
    auto e2 = boost::add_edge(v, u, er, g_).first;
    g_[e1].reverse_edge = e2;
    g_[e2].reverse_edge = e1;
}


double GraphModule::max_flow(const int source, const int sink)
{
    auto idx = get(&Atom::id, g_);
    auto cap    = get(&HydrogenBond::energy, g_);
    auto rescap = get(&HydrogenBond::residual_energy, g_);
    auto rev = get(&HydrogenBond::reverse_edge, g_);
	  
    double flow = boost::boykov_kolmogorov_max_flow(g_, cap, rescap, rev, idx,
						    source, sink);
    return flow;
}

void GraphModule::clear()
{
    g_.clear();
}

