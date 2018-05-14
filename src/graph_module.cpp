
#include "graph_module.hpp"
#include <boost/make_shared.hpp>
#include <boost/range/iterator_range.hpp>


GraphModule::GraphModule(int n)
{
    this->g_ = std::make_unique<graph_t>(graph_t(n));
}

graph_t GraphModule::get_graph() const
{
    return *this->g_;
}

bool GraphModule::add_bidirectional_edge(int i, int j, const EdgeProperties& properties)
{
    graph_t::edge_descriptor edge;
    bool inserted;
    
    boost::tie(edge, inserted) = add_edge(i, j, properties, *this->g_);
    boost::tie(edge, inserted) = add_edge(j, i, properties, *this->g_);

    return inserted;
}

// float GraphModule::flow(int i, int j)
// {
//     graph_t::vertex_descriptor s = vertex(i, *this->g_);
//     graph_t::vertex_descriptor t = vertex(j, *this->g_);
//     // return boost::boykov_kolmogorov_max_flow(*this->g_,
//     // 					     get(&EdgeProperties::capacity, *this->g_),
//     // 					     get(&EdgeProperties::residual_capacity, *this->g_),
//     // 					     get(&EdgeProperties::reverse, *this->g_),
//     // 					     get(&VertexProperties::id, *this->g_), s, t);
//     return boost::boykov_kolmogorov_max_flow(*this->g_,
// 					     capacity_map(get(&EdgeProperties::capacity, *this->g_),
// 					     get(&EdgeProperties::residual_capacity, *this->g_),
// 					     get(&EdgeProperties::reverse, *this->g_),
// 					     get(boost::vertex_index, *this->g_), s, t);
//     return 2.0;
// }

/*
std::vector<ComponentGraph> connected_components_subgraphs(const Graph &g)
{
    vertex_component_map mapping =
	boost::make_shared<std::vector<unsigned long>>(num_vertices(g));
    size_t num = boost::connected_components(g, mapping->data());

    std::vector<ComponentGraph> component_graphs;

    for (size_t i = 0; i < num; i++)
	component_graphs.emplace_back(g,
				      [mapping,i,&g](Graph::edge_descriptor e) {
					  return mapping->at(source(e,g))==i
					      || mapping->at(target(e,g))==i;}, 
				      [mapping,i](Graph::vertex_descriptor v) {
					  return mapping->at(v)==i;});
    
    // for (auto &comp : component_graphs) {
    // 	std::cout << boost::distance(boost::vertices(comp)) << " ";
    // }
    
    // std::cout << "\n";
    
    return component_graphs;
}

Subgraph buildSubgraph(const std::vector<int> &vertices, Subgraph& g) {
    Subgraph& s = g.create_subgraph();
    for (auto vertex : vertices) {
	add_vertex(vertex, s);
    }
    return s;
}
*/
