
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from itertools import combinations


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
    print('Read {} Graph nodes'.format(len(nodes)))
    return nodes

def analyze_subgraph(g, n):
    sub = nx.subgraph(g, \
                      nx.node_connected_component(g, n))
    p = nx.current_flow_betweenness_centrality(sub,normalized=False,weight='weight')
    avg = nx.global_efficiency(sub)
    return sub, avg, len(sub), p[n]

if __name__ == "__main__":
    
    print('NetworkX Version : {}'.format(nx.__version__))
    
    frames = None
    sizes = []
    current_flow_betweenness_centrality = []
    global_efficiency = []
    subgraphs = []
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
                        subgraphs.append(sub)
                         
                        p = nx.current_flow_betweenness_centrality(sub,\
                                                                   normalized=False,\
                                                                   weight='weight')
                        sizes.append(len(sub))
                        current_flow_betweenness_centrality.append(p[n])
                        global_efficiency.append(nx.global_efficiency(sub))
                G = nx.Graph() # Reset
                print('{} frames read found {} networks... '.format(current, len(subgraphs)), end='\r')
                continue
            # if len(subgraphs) > 100:
            #     break
            if current == -1:
                break
            
            line = line.split()
            G.add_edge(int(line[0]), int(line[1]), \
                       weight=abs(float(line[2])), \
                       length=float(line[3]), angle=float(line[4]))
            
    print('\nfound {} networks'.format(len(subgraphs)))
    
    ids = range(len(subgraphs))
    count = 0
    iso_pairs = []
    combs = len(combinations(ids, 2))
    for i, j in combinations(ids, 2):
        count += 1
        if (count % combs/10) == 0:
            print('{} done'.format(int(count/combs*100)), end='\r')
        if nx.is_isomorphic(subgraphs[i], subgraphs[j]):
            iso_pairs.append((i, j))
    
    iso_graph = nx.from_edgelist(iso_pairs)
    
    # for c in nx.connected_components(iso_graph):
    #    print(c)
    iso_cluster = []
    iso_sizes = []
    for c in nx.connected_components(iso_graph):
        avg = 0
        for i in c:
            avg += sizes[i]
        avg /= len(c)
        iso_cluster.append(len(c))
        iso_sizes.append(avg)
    
    plt.scatter(iso_sizes, iso_cluster)
    plt.xlabel('sizes')
    plt.ylabel('cluster')
    # nx.draw(iso_graph, with_labels=True)
    
    # id_max = np.argmax(sizes)
    # id_min = np.argmin(sizes)
    # print_graph(subgraphs[id_max], 'test.png')

    # fig1, ax1 = plt.subplots()
    # h1 = ax1.scatter(sizes, global_efficiency, c=current_flow_betweenness_centrality, \
    #                  s=100, edgecolor='') 
    # ax1.set_xlabel('sizes')
    # ax1.set_ylabel('global efficiency')

    # fig2, ax2 = plt.subplots()
    # h2 = ax2.scatter(sizes, current_flow_betweenness_centrality, s=100, edgecolor='') 
    # ax2.set_xlabel('sizes')
    # ax2.set_ylabel('centrality')
    
    # fig3, ax3 = plt.subplots()
    # h3 = ax3.scatter(global_efficiency, current_flow_betweenness_centrality, \
    #                  s=100, edgecolor='') 
    # ax3.set_xlabel('avgs')
    # ax3.set_ylabel('centrality')
    
    plt.show()
    
    
