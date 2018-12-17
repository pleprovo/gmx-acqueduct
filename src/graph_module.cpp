
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


void GraphModule::add_edge1(const int u, const int v, const HydrogenBond &e)
{
    auto e1 = boost::add_edge(u, v, e, g_).first;
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

    // boykov_kolmogorov_max_flow(Graph& g,
    //    CapacityEdgeMap cap,
    //    ResidualCapacityEdgeMap res_cap,
    //    ReverseEdgeMap rev_map,
    //    PredecessorMap pre_map,
    //    ColorMap color,
    //    DistanceMap dist,
    //    IndexMap idx,
    //    typename graph_traits <Graph>::vertex_descriptor src,
    //    typename graph_traits <Graph >::vertex_descriptor sink)
    
}

std::vector<int> GraphModule::dfs(const int source)
{
    std::vector<int> predecessor_map;

    auto indexmap = boost::get(&Atom::id, g_);
    auto colormap = boost::make_vector_property_map<boost::default_color_type>(indexmap);
    std::vector<int> p(num_vertices(g_));
    std::vector<int> d(num_vertices(g_));
    
    boost::depth_first_search(g_, vis_, colormap, source);
     // boost::breadth_first_search(g_, source, 
     // 				boost::make_bfs_visitor(
     // 				    std::make_pair(boost::record_distances(d, boost::on_tree_edge()), boost::record_predecessors(p.begin(), boost::on_tree_edge()))));
    // std::vector<int> vctr = vis_.GetVector();

    // for(auto id : vctr)
    //     std::cout << id << " ";
    return vis_.GetVector();
}

void GraphModule::clear()
{
    g_.clear();
}

