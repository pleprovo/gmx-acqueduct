import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import cm
from itertools import combinations
from math import factorial, ceil, log
import multiprocessing
from joblib import Parallel, delayed
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import concurrent.futures


def filter_edge_by_weight(G, threshold=-1):
    """
    Remove all edge below threshold in the graph G.
    """
    g2 = G.copy()
    for edge in g2.edges():
        if G[edge[0]][edge[1]]['weight'] <= threshold:
            G.remove_edge(edge[0], edge[1])
    return G

def remove_node_degree(G, degree=0):
    """
    Remove all nodes without edges in the graph G.
    """
    degrees = G.degree()
    g2 = G.copy()
    for node in g2.nodes():
        if degrees[node] == degree:
            G.remove_node(node)
    return G

def print_graph(g, filename):
    """
    Print a fancy graph to file.
    """
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
    """
    Read Node information from file.
    """
    nodes = dict()
    resis = dict()
    next(infile)
    for line in infile:
        line = line.split()
        nodes[int(line[0])] = "{} {} {}".format(line[2].strip(), \
                                                line[3].strip(), \
                                                line[4].strip())
        resis['{} {}'.format(line[2].strip(), line[3].strip())] = int(line[0])
    print('Read {} Graph nodes.'.format(len(nodes)))
    return nodes, resis

def isomorph_pair_graph(frames, method=nx.faster_could_be_isomorphic):
    """
    Find Possible Isomorph Pair Between each network in frames.
    """
    ids = range(frames.shape[0])
    count = 0
    num_comb = factorial(len(ids))/(factorial(2)*factorial(len(ids)-2))
    iso_pairs = []
    print('0% combinations done ...', end='\r')
    for i, j in combinations(ids, 2):
        count += 1
        if (count % int(num_comb/100)) == 0:
            print('{}% combinations done ...'.format(ceil(count/num_comb*100)), end='\r')
        if method(frames['graph'][i], frames['graph'][j]):
            iso_pairs.append((i, j))
    
    print('{} combinations done.'.format(count))
    return nx.from_edgelist(iso_pairs)

def analyze_isomorph_components(frames, isomorph_graph):
    """
    Analysis average characteristic over isomorph clusters (connected component \
    in isomorph_graph)
    """
    count = 0
    results = []
    for i, c in enumerate(nx.connected_components(isomorph_graph)):
        print('{} components done ...'.format(count), end='\r')
        num_edge = 0
        efficiency = 0
        centrality = 0
        for elem in c:
            num_edge += frames['num_edges'][elem]
            centrality += frames['centrality'][elem]
            efficiency += frames['efficiency'][elem]

        new_row = {'id' : i, \
                   'size' : len(c), \
                   'avg_size' : num_edge, \
                   'centrality' : centrality, \
                   'efficiency' :  efficiency}
        results.append(new_row)
    return pd.DataFrame(results)

def pca_frames(frames, features=['num_nodes', 'num_edges', \
                                 'centrality', 'efficiency', \
                                 'clustering'], normalized=True, num_components=2):
    x = frames.loc[:, features].values
    
    if normalized:
        x = StandardScaler().fit_transform(x)

    pca = PCA(n_components=num_components)

    principalComponents = pca.fit_transform(x)
    
    return pd.DataFrame(data = principalComponents
                               , columns = ['principal component 1', 'principal component 2'])

def main():

    nodes = None
    residues = None
    with open('nodes-list.xvg.dat', 'r') as infile:
        nodes, residues = read_node_list(infile)

    print(residues['SG 374'], residues['SG 763'], residues['SG 1152'], residues['SG 1541'])
        
    data = []
    with open('edges-frames.xvg.dat', 'r') as infile:
        next(infile)
        current = 0
        G = nx.Graph()
        for line in infile:
            if line[0] == '#':
                current += 1
                G = filter_edge_by_weight(G, threshold=0.0)
                G = remove_node_degree(G)
                for i, n in enumerate(['SG 374', 'SG 763', 'SG 1152', 'SG 1541']):
                    id = residues[n]
                    if id in G:
                        sub = nx.subgraph(G, nx.node_connected_component(G, id))
                        nx.set_node_attributes(sub, name='label', values=nodes)
                        new_network = {'num_frames' : current, \
                                       'chain' : i, \
                                       'num_edges' : nx.number_of_edges(sub), \
                                       'num_nodes' : nx.number_of_nodes(sub), \
                                       'centrality' : nx.current_flow_betweenness_centrality(sub, normalized=False, weight='weight')[id], \
                                       'efficiency' : nx.global_efficiency(sub), \
                                       'clustering' : nx.average_clustering(sub), \
                                       'avg_connectivity' : nx.average_node_connectivity(sub), \
                                       'avg_shortest_path' : nx.average_shortest_path_length(sub), \
                                       'graph' : sub }
                        data.append(new_network)
                        G = nx.Graph()
                    print('{} frames read and {} networks found ...'.format(current, len(data)), end='\r')
                continue
            if len(data) >= 100:
                break
            if current == -1:
                break
            
            line = line.split()
            G.add_edge(int(line[0]), int(line[1]), \
                       weight=abs(float(line[2])), \
                       length=float(line[3]), angle=float(line[4]))
        print('{} frames read.'.format(current))

    frames = pd.DataFrame(data)
    
    frames.num_edges = frames.num_edges.astype(np.int32)
    frames.num_nodes = frames.num_nodes.astype(np.int32)
    frames.num_frames = frames.num_frames.astype(np.int32)
    frames.centrality = frames.centrality.astype(np.float32)
    frames.clustering = frames.clustering.astype(np.float32)
    frames.efficiency = frames.efficiency.astype(np.float32)
    frames.avg_connectivity = frames.avg_connectivity.astype(np.float32)
    frames.avg_shortest_path = frames.avg_shortest_path.astype(np.float32)
    
    del data[:]
    print('{} network found.'.format(frames.shape[0]))
    
    ## ISOMORPHISM ##
    isomorph_graph = isomorph_pair_graph(frames)

    avg_isomorph_components = analyze_isomorph_components(frames, isomorph_graph)

    ## PCA SCIKIT ##
    features = ['num_edges', 'num_nodes', 'centrality', 'clustering', 'efficiency', 'avg_connectivity', 'avg_shortest_path']

    finalDf = pca_frames(frames)

    figs = []
    axes = []
    handle = []
    
    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1) 
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_title('2 Component PCA -', fontsize = 20)

    h = ax.scatter(finalDf.loc[:, 'principal component 1'], \
                   finalDf.loc[:, 'principal component 2'], \
                   s = 50)
    ax.grid()
    figs.append(fig)
    axes.append(ax)
    handle.append(h)

    fig, ax = plt.subplots()
    h = ax.hist(frames.chain, bins=4)
    ax.set_xlabel('Chains')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Number of Networks per Chain')
    figs.append(fig)
    axes.append(ax)
    handle.append(h)
    
    fig, ax = plt.subplots()
    h = ax.hist(frames.num_edges, bins=50)
    ax.set_xlabel('Sizes' )
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Networks Sizes')
    figs.append(fig)
    axes.append(ax)
    handle.append(h)

    # fig, ax = plt.subplots()
    # h = ax.hist(frames.centrality, bins=50)
    # ax.set_xlabel('Centrality')
    # ax.set_ylabel('Count')
    # ax.set_title('Distribution of Networks Centrality')
    # figs.append(fig)
    # axes.append(ax)
    # handle.append(h)

    # fig, ax = plt.subplots()
    # h = ax.hist(frames.efficiency, bins=50)
    # ax.set_xlabel('Efficiency')
    # ax.set_ylabel('Count')
    # ax.set_title('Distribution of Networks Efficiency')
    # figs.append(fig)
    # axes.append(ax)
    # handle.append(h)

    # fig, ax = plt.subplots()
    # h = ax.hist(frames.clustering, bins=50)
    # ax.set_xlabel('Average Clustering ')
    # ax.set_ylabel('Count')
    # ax.set_title('Distribution of Network Average Clustering')
    # figs.append(fig)
    # axes.append(ax)
    # handle.append(h)

    # fig, ax = plt.subplots()
    # h = ax.hist(frames.avg_connectivity, bins=50)
    # ax.set_xlabel('Average Node Connectivity')
    # ax.set_ylabel('Count')
    # ax.set_title('Distribution of Networks Average Node Connectivity')
    # figs.append(fig)
    # axes.append(ax)
    # handle.append(h)

    # fig, ax = plt.subplots()
    # h = ax.hist(frames.avg_shortest_path, bins=50)
    # ax.set_xlabel('Average Shortest Path Length')
    # ax.set_ylabel('Count')
    # ax.set_title('Distribution of Average Sortest Path Length')
    # figs.append(fig)
    # axes.append(ax)
    # handle.append(h)
    
    

    # fig, ax = plt.subplots()
    # h = ax.scatter(frames.num_nodes, frames.efficiency, \
    #                c=frames.centrality, \
    #                s=100) 
    # ax.set_xlabel('sizes')
    # ax.set_ylabel('global efficiency')
    # ax.set_title('Centrality on sizes')
    # figs.append(fig)
    # axes.append(ax)
    # handle.append(h)
    
    # fig, ax = plt.subplots()
    # h = ax.scatter(frames.num_nodes, frames.centrality, \
    #                  c=frames.efficiency, \
    #                  s=100) 
    # ax.set_xlabel('sizes')
    # ax.set_ylabel('centrality')
    # ax.set_title('Efficiency on sizes')
    # figs.append(fig)
    # axes.append(ax)
    # handle.append(h)
    
    # fig, ax = plt.subplots()
    # h = ax.scatter(frames.efficiency, frames.centrality, \
    #                c=frames.num_nodes, s=100) 
    # ax.set_xlabel('avgs')
    # ax.set_ylabel('centrality')
    # figs.append(fig)
    # axes.append(ax)
    # handle.append(h)

    fig, ax = plt.subplots()
    h = ax.scatter(avg_isomorph_components.centrality, avg_isomorph_components.efficiency, \
                   s=log(avg_isomorph_components.size+1), \
                   alpha=0.5)
    plt.colorbar(h)
    ax.set_xlabel('Centrality')
    ax.set_ylabel('Efficiency')
    ax.set_title('Centratility on Sizes')
    figs.append(fig)
    axes.append(ax)
    handle.append(h)
    
    # sizes_networks = np.array(sizes_networks).T
    # fig, ax = plt.subplots()
    # h1 = ax.plot(sizes_networks[0])
    # h2 = ax.plot(sizes_networks[1])
    # h3 = ax.plot(sizes_networks[2])
    # h4 = ax.plot(sizes_networks[3])
    # ax.set_xlabel('frames')
    # ax.set_ylabel('Size of networks detected')
    # figs.append(fig)
    # axes.append(ax)
    # handle.append(h1)
    # handle.append(h2)
    # handle.append(h3)
    # handle.append(h4)
    
    plt.show()


    
if __name__ == "__main__":
    
    print('NetworkX Version : {}'.format(nx.__version__))
    main()
    
