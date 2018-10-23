# gmx-waternetwork
Gromacs Analysis tools for Water Hydrogen Bond Networks consisting in two programs for analysing buried waters:
1. Alpha Shape tool to compute the popuplation of water buried within a surface define by th user
2. Waternetwork tool that use triangulation and an hydrogen bond potential to define a hydrogen 
bond network. Graph characterisic can be compute on the fly (shortest path, maximum flow, cluster, isomorphisme)
## Dependencies
1. boost (the most recent vesion)
2. GROMACS 5.1.2
3. CGAL 4.X
## Build instructions
1. mkdir build
2. cd build/
3. cmake ..
4. make -k
## Usage


