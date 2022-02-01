import os
from network_vis import *
from pymol import cmd

def load_network(points, edges):
    with open(points, 'r') as infile:
        for line in infile.readlines():
            field = line.strip().split()
            xyz = map(float, field)
            xyz = [10*x for x in xyz]
            cmd.pseudoatom("test", elem='O', pos=xyz)

    with open(edges, 'r') as infile:
        for line in infile.readlines():
            field = line.strip().split()
            if float(field[-1]) > 0.1:
                plop = cmd.bond('index %s' % field[0], 'index %s' % field[1])
                cmd.set_bond('line_color', float(field[-1]), 'i. 1')

    cmd.show("sphere", "test")
    cmd.set('line_width', 5.0)
    return 0
    
cmd.extend("load_network", load_network)
