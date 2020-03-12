
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import cm
from itertools import combinations
from math import factorial, ceil
import multiprocessing
from joblib import Parallel, delayed

def filter_edge_by_weight(G, threshold=-1):
    g2 = G.copy()
    for edge in g2.edges():
        if G[edge[0]][edge[1]]['weight'] <= threshold:
            G.remove_edge(edge[0], edge[1])
    return G

def remove_node_degree(G, degree=0):
    degrees = G.degree()
    g2 = G.copy()
    for node in g2.nodes():
        if degrees[node] == degree:
            G.remove_node(node)
    return G

def print_graph(g, filename):
    node_colors = nx.get_node_attributes(g,'color').values()
    edge_weights = nx.get_edge_attributes(g,'weight').values()
    emin = min(edge_weights)
    emax = max(edge_weights)
    
    fig, ax = plt.subplots()
    pos = nx.spring_layout(g, weight='length')
    
    nx.draw_networkx_nodes(g, pos, nodelist=g.nodes(), node_size=50)
        
    nx.draw_networkx_edges(g, pos, edgelist=g.edges(), \
                           width=5, edge_cmap=cm.viridis, \
                           edge_vmax=emax, edge_vmin=emin, edge_color=edge_weights)
    
    nx.draw_networkx_labels(g,pos, nx.get_node_attributes(g,'label'))
    
    plt.savefig(filename)
    # plt.cla()
    # plt.clf()
        
def read_node_list(infile):
    nodes = dict()
    next(infile)
    for line in infile:
        line = line.split()
        nodes[int(line[0])] = "{} {} {}".format(line[2].strip(), \
                                                line[3].strip(), \
                                                line[4].strip())
    print('Read {} Graph nodes.'.format(len(nodes)))
    return nodes

if __name__ == "__main__":
    
    print('NetworkX Version : {}'.format(nx.__version__))
    
    frames = None
    subgraphs = []
    sizes = []
    current_flow_betweenness_centrality = []
    global_efficiency = []
    nodes = None
    
    with open('node_list.txt', 'r') as infile:
        nodes = read_node_list(infile)

    with open('edge_list.txt', 'r') as infile:
        next(infile)
        current = 0
        G = nx.Graph()
        for line in infile:
            if line[0] == '#':
                current += 1
                # frames[current] = nx.Graph()
                G = filter_edge_by_weight(G, threshold=0.0)
                G = remove_node_degree(G)
                for n in [959, 1950, 2941, 3932]:
                    if n in G:
                        sub = nx.subgraph(G, nx.node_connected_component(G, n))
                        nx.set_node_attributes(sub, name='label', values=nodes)
                        p = nx.current_flow_betweenness_centrality(sub, \
                                                   normalized=False, weight='weight')
                        sizes.append(len(sub))
                        current_flow_betweenness_centrality.append(p[n])
                        global_efficiency.append(nx.global_efficiency(sub))  
                        subgraphs.append(sub)
                G = nx.Graph() # Reset
                print('{} frames read ...'.format(current, len(subgraphs)), end='\r')
                continue
            # if len(subgraphs) >= 250:
            #     break
            if current == 1000:
                break
            
            line = line.split()
            G.add_edge(int(line[0]), int(line[1]), \
                       weight=abs(float(line[2])), \
                       length=float(line[3]), angle=float(line[4]))
        print('{} frames read.'.format(current))
        
    print('{} network found.'.format(len(subgraphs)))
    
    ids = range(len(subgraphs))
    count = 0
    num_comb = factorial(len(ids))/(factorial(2)*factorial(len(ids)-2))
    
    iso_pairs = []
    print('0% combinations done ...', end='\r')
    for i, j in combinations(ids, 2):
        count += 1
        if (count % int(num_comb/100)) == 0:
            print('{}% combinations done ...'.format(ceil(count/num_comb*100)), end='\r')
        if nx.is_isomorphic(subgraphs[i], subgraphs[i]):
            iso_pairs.append((i, j))
            
    print('{} combinations done.'.format(count))
    iso_graph = nx.from_edgelist(iso_pairs)

    cluster_sizes = []
    iso_efficiency = []
    iso_sizes = []
    iso_centrality = []
    count = 0
    for c in nx.connected_components(iso_graph):
        count += 1
        print('{} components done ...'.format(count), end='\r')
        cluster_sizes.append(len(c))
        size = 0
        efficiency = 0
        centrality = 0
        for elem in c:
            size += sizes[elem]
            centrality += current_flow_betweenness_centrality[elem]
            efficiency += global_efficiency[elem]
        iso_sizes.append(size/cluster_sizes[-1])
        iso_efficiency.append(efficiency/cluster_sizes[-1])        
        iso_centrality.append(centrality/cluster_sizes[-1])

    figs = []
    axes = []
    handle = []

    fig, ax = plt.subplots()
    h = ax.hist(sizes, bins=50)
    ax.set_xlabel('Sizes')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Networks Sizes')
    figs.append(fig)
    axes.append(ax)
    handle.append(h)

    fig, ax = plt.subplots()
    h = ax.hist(current_flow_betweenness_centrality, bins=50)
    ax.set_xlabel('Centrality')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Networks Centrality')
    figs.append(fig)
    axes.append(ax)
    handle.append(h)

    fig, ax = plt.subplots()
    h = ax.hist(global_efficiency, bins=50)
    ax.set_xlabel('Efficiency')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Networks Efficiency')
    figs.append(fig)
    axes.append(ax)
    handle.append(h)
    
    fig, ax = plt.subplots()
    h = ax.scatter(iso_centrality, iso_efficiency, \
                   c=iso_sizes, \
                   s=cluster_sizes, \
                   alpha=0.5)
    plt.colorbar(h)
    ax.set_xlabel('cluster sizes')
    ax.set_ylabel('Average Network Size')
    ax.set_title('Centratility on Sizes')
    figs.append(fig)
    axes.append(ax)
    handle.append(h)

    # fig, ax = plt.subplots()
    # h = ax.scatter(cluster_sizes, iso_sizes, \
    #                  c=iso_efficiency, \
    #                  s=100)
    # plt.colorbar(h)
    # ax.set_xlabel('cluster sizes')
    # ax.set_ylabel('Average Network Size')
    # ax.set_title('Efficiency on Sizes')
    # figs.append(fig)
    # axes.append(ax)
    # handle.append(h)
    
    # fig, ax = plt.subplots()
    # h = ax.scatter(iso_centrality, iso_efficiency, \
    #                  c=iso_sizes, \
    #                  s=100)
    # plt.colorbar(h)
    # ax.set_xlabel('Centrality')
    # ax.set_ylabel('Efficiency')
    # ax.set_title('Sizes on measures')
    # figs.append(fig)
    # axes.append(ax)
    # handle.append(h)

    # fig, ax = plt.subplots()
    # h = ax.scatter(sizes, global_efficiency, \
    #                  c=current_flow_betweenness_centrality, \
    #                  s=100) 
    # ax.set_xlabel('sizes')
    # ax.set_ylabel('global efficiency')
    # ax.set_title('Centrality on sizes')
    # figs.append(fig)
    # axes.append(ax)
    # handle.append(h)
    
    # fig, ax = plt.subplots()
    # h = ax.scatter(sizes, current_flow_betweenness_centrality, \
    #                  c=global_efficiency, \
    #                  s=100) 
    # ax.set_xlabel('sizes')
    # ax.set_ylabel('centrality')
    # ax.set_title('Efficiency on sizes')
    # figs.append(fig)
    # axes.append(ax)
    # handle.append(h)
    
    # fig, ax = plt.subplots()
    # h = ax.scatter(global_efficiency, current_flow_betweenness_centrality, \
    #                  c=sizes, s=100) 
    # ax.set_xlabel('avgs')
    # ax.set_ylabel('centrality')
    # figs.append(fig)
    # axes.append(ax)
    # handle.append(h)
    
    plt.show()
    
    
