
#include "graph_module.hpp"
#include <boost/make_shared.hpp>
#include <boost/range/iterator_range.hpp>
#include <iostream>


GraphModule::GraphModule () {

}

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

float averageDegree(const Graph &g) {
    vertex_iterator vi, vi_end, next;
    boost::tie(vi, vi_end) = boost::vertices(g);
    float average = 0.0;
    for (next = vi; vi != vi_end; vi = next) {
	++next;
        average += boost::degree(*vi, g);
    }
    return average/boost::num_vertices(g);

}

int writeGraphEdges(const Graph& g, std::stringstream& oss) {
    oss << "Edges " << "\n";
    // Iterate through the edges and print them out
    Graph::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
	oss << boost::source(*ei,g) << " " << boost::target(*ei, g) << "\n";
    }
    return 0;
}
