
#include "graph_module.hpp"


void print_predecessor_path(Graph &g, Traits::vertex_descriptor v)
{
    using path_t = std::vector<Graph::edge_descriptor>;
    path_t path;    
    for(Graph::vertex_descriptor u = g[v].predecessor; u != v; v=u, u=g[v].predecessor)
    {
    	std::pair<Graph::edge_descriptor, bool> edge_pair = boost::edge(u,v,g);
    	path.push_back( edge_pair.first );
    }
        
    std::cout << "Shortest Path from v1 to v6:" << std::endl;
    for(path_t::reverse_iterator riter = path.rbegin(); riter != path.rend(); ++riter)
    {
        Graph::vertex_descriptor u_tmp = boost::source(*riter, g);
        Graph::vertex_descriptor v_tmp = boost::target(*riter, g);
        Graph::edge_descriptor e_tmp = boost::edge(u_tmp, v_tmp, g).first;
	
    	std::cout << "  " << g[u_tmp].id << " -> " << g[v_tmp].id << "    (weight: " << g[e_tmp].length << ")" << std::endl;
    }
}

double do_max_flow(Graph &g, Graph::vertex_descriptor &source, Graph::vertex_descriptor &sink)
{
    auto idx = get(&Atom::id, g);
    auto cap    = get(&HydrogenBond::energy, g);
    auto rescap = get(&HydrogenBond::residual_energy, g);
    auto rev = get(&HydrogenBond::reverse_edge, g);

    double flow = boost::boykov_kolmogorov_max_flow(g, cap, rescap, rev, idx, source, sink);
    return flow;
}

GraphModule::GraphModule()
{
    g_ = Graph();
}

GraphModule::GraphModule(const int n)
{
    g_ = Graph(n);
}

void GraphModule::add_vertex(const Atom &v)
{
    boost::add_vertex(v, g_);
}


void GraphModule::add_edge(const int u, const int v, const HydrogenBond &e)
{
    auto e1 = boost::add_edge(u, v, e, g_).first;
    auto e2 = boost::add_edge(v, u, e, g_).first;
    g_[e1].reverse_edge = e2;
    g_[e2].reverse_edge = e1;
}


double GraphModule::max_flow(const int source, const int sink)
{
    auto idx = get(&Atom::id, g_);
    auto cap    = get(&HydrogenBond::energy, g_);
    auto rescap = get(&HydrogenBond::residual_energy, g_);
    auto rev = get(&HydrogenBond::reverse_edge, g_);

    double flow = boost::boykov_kolmogorov_max_flow(g_, cap, rescap, rev, idx, source, sink);
    return flow;
}


void GraphModule::clear()
{
    g_.clear();
}

// void addBidirectionalEdge(Graph &graph,
// 			  const unsigned int &source,
// 			  const unsigned int &target,
// 			  const float &weight,
//                           std::vector<Traits::edge_descriptor>& reverseEdges,
// 			  std::vector<float>& capacity)
// {
//     // Add edges between grid vertices. We have to create the edge and the reverse edge,
//     // then add the reverseEdge as the corresponding reverse edge to 'edge', and then add 'edge'
//     // as the corresponding reverse edge to 'reverseEdge'
    
//     //int nextEdgeId = num_edges(graph);
//     int nextEdgeId = capacity.size();
//     Traits::edge_descriptor edge;
//     bool inserted;

//     boost::tie(edge,inserted) = add_edge(source, target, nextEdgeId, graph);
//     if(!inserted) {
//         std::cerr << "Not inserted!" << std::endl;
//     }
//     Traits::edge_descriptor reverseEdge = add_edge(target, source, nextEdgeId + 1, graph).first;
//     reverseEdges.push_back(reverseEdge);
//     reverseEdges.push_back(edge);
//     capacity.push_back(weight);

//     // Not sure what to do about reverse edge weights
//     capacity.push_back(weight);
//     // capacity.push_back(0);
// }

// GraphAnalysisModule::GraphAnalysisModule()
// {
//     g_ = Network();
// }

// void GraphAnalysisModule::add_atom(int id, int resid, std::string name, bool type)
// {
//     add_vertex(Atoms{id, resid, name, type}, g_);
// }

// void GraphAnalysisModule::add_bidirectional_edge(int v1, int v2, HydrogenBonds &bond)
// {
//     auto a = boost::add_edge(v1, v2, bond, g_).first;
//     auto b = boost::add_edge(v2, v1, bond, g_).first;
//     reverse_edge_[a] = b;
//     reverse_edge_[b] = a;
// }

// float GraphAnalysisModule::flow(int source, int target)
// {
//     struct VertexEx {
//     	boost::default_color_type color;
//         double distance;
//         Graph::edge_descriptor pred;
//     };

//     auto idx = get(boost::vertex_index, g_);
//     auto vex = boost::make_vector_property_map<VertexEx>(idx);
    
//     auto pred = boost::make_transform_value_property_map(std::mem_fn(&VertexEx::pred), vex);
//     auto color = boost::make_transform_value_property_map(std::mem_fn(&VertexEx::color), vex);
//     auto dist = boost::make_transform_value_property_map(std::mem_fn(&VertexEx::distance), vex);

//     auto weight = get(&HydrogenBonds::length, g_);
//     auto cap    = get(&HydrogenBonds::capacity, g_);
//     auto rescap = get(&HydrogenBonds::residual_capacity, g_);

//     double flow = boost::boykov_kolmogorov_max_flow(g_, cap, rescap, reverse_edge_,
// 						    pred, color, dist, idx, source, target);

//     return flow;
// }

