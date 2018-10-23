
#include <iostream>

#include "graph_module.hpp"


int main (int argc, char *argv[])
{
    std::cout << " --- TESTS --- " << std::endl;

    Graph g;
    GraphModule gm;
    Atom a1{0};
    
    gm.add_vertex(Atom{0});
    gm.add_vertex(Atom{1});
    gm.add_vertex(Atom{2});
    gm.add_vertex(Atom{3});
    gm.add_vertex(Atom{4});
    gm.add_vertex(Atom{5});

    gm.add_edge(0, 1, HydrogenBond{3.5, 0.0, 4.0});
    gm.add_edge(1, 2, HydrogenBond{3.5, 0.0, 4.0});
    gm.add_edge(1, 3, HydrogenBond{3.5, 0.0, 4.0});
    gm.add_edge(3, 4, HydrogenBond{3.5, 0.0, 4.0});
    gm.add_edge(2, 4, HydrogenBond{3.5, 0.0, 4.0});
    gm.add_edge(4, 5, HydrogenBond{3.5, 0.0, 4.0});

    // std::cout <<  " Test Graph : " << std::endl;
    // boost::print_graph(g, boost::get(&Atom::resid, g));
    // std::cout << std::endl;

    // boost::dijkstra_shortest_paths(g, v0, boost::weight_map(get(&HydrogenBond::length,g)).distance_map(get(&Atom::distance, g)).predecessor_map(get(&Atom::predecessor, g)));
    
    // print_predecessor_path(g, v5);
    

    double flow = gm.max_flow(0, 5);
    std::cout << "Flow: " << flow << "\n";
    gm.clear();

    // Testing FFTW
    // int N = data.size();
    // double in[N];
    // fftw_complex *out;
    // fftw_plan fft, ifft;

    // out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
    // fft = fftw_plan_dft_r2c_1d(N, in, out, FFTW_FORWARD);
    
    // fftw_execute(my_plan);    
    // fftw_destroy_plan(my_plan);
    
    // fftw_free(in);    
    // fftw_free(out);
    
    return 0;
}
