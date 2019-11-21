# gmx-waternetwork

Gromacs Analysis tools for Water Hydrogen Bond Networks consisting in two programs for analysing buried waters:

* Alpha Shape tool to compute the popuplation of water buried within a surface define by th user
* Waternetwork tool that use triangulation and an hydrogen bond potential to define a hydrogen 
bond network. Graph characterisic can be compute on the fly (shortest path, maximum flow, cluster, isomorphisme)

## Dependencies
1. [boost](https://www.boost.org/doc/libs/1_70_0/libs/graph/doc/) (the most recent vesion)
2. [GROMACS 5.1.2](http://manual.gromacs.org/documentation/5.1.2/install-guide/index.html)
3. [CGAL 4.X](https://www.cgal.org/)

## Build instructions
```
mkdir build
cd build/
cmake ..
make -k
```

## Usage
```
gmx waternetwork -f *.xtc -s * -o all.xvg -n *.ndx
```

