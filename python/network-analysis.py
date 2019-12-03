
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import cm

print('NetworkX Version : {}'.format(nx.__version__))

def get_minmax_edges(G, keyword):
    edge_data = nx.get_edge_attributes(G, keyword).values()
    emax = max(edge_data)
    emin = min(edge_data)
    print('max : {}, min : {}'.format(emax, emin))
    return emin, emax, edge_data

def filter_edge_by_weight(G, threshold=-1):
    g2 = G.copy()
    for edge in g2.edges():
        if G[edge[0]][edge[1]]['weight'] < threshold:
            G.remove_edge(edge[0], edge[1])
    return G

def remove_node_degree(G, degree=0):
    degrees = G.degree()
    g2 = G.copy()
    for node in g2.nodes():
        if degrees[node] == degree:
            G.remove_node(node)
    return G

def print_graph(D, filename):
    emin, emax, edge_weigths = get_minmax_edges(D, 'weight')

    pos = nx.kamada_kawai_layout(D)
    nx.draw_networkx_nodes(D, pos, nodelist=D.nodes(), \
                           node_size=200, node_Golor=nx.get_node_attributes(D,'color').values())
        
    nx.draw_networkx_edges(D, pos, edgelist=D.edges(), \
                           width=5, edge_cmap=cm.hot, \
                           edge_vmax=emax, edge_vmin=emin, edge_Golor=nx.get_edge_attributes(D,'weight').values())
    nx.draw_networkx_labels(D,pos, nx.get_node_attributes(D,'label'))
    plt.savefig(filename)
    plt.cla()
    plt.clf()

if __name__ == "__main__":
    G = nx.Graph()
    with open('node_list.txt', 'r') as infile:
        next(infile)
        for line in infile:
            line = line.split()
            color = ()
            if "WAT" in line:
                color = 'r'
            elif "SCY" in line:
                color = 'g'
            else:
                color = 'b'
            G.add_node(int(line[0]), label="{} {} {}".format(line[2].strip(), \
                                                             line[3].strip(), \
                                                             line[4].strip()), \
                       color=color)    

    with open('edge_list.txt', 'r') as infile:
        next(infile)
        next(infile)
        count = 0
        for line in infile:
            count += 1
            if line[0] == '#':
                break
            # if count % 1000 == 0:
            #    print('{} edges'.format(count))
            line = line.split()
            G.add_edge(int(line[0]), int(line[1]), weight=abs(float(line[2])), \
                       length=float(line[3]), angle=float(line[4]))

    print('Num of nodes : {}, Num of edges : {}'.format(G.number_of_nodes(), G.size()))
    G = filter_edge_by_weight(G, threshold=0.5)
    G = remove_node_degree(G)  
    print('Num of nodes : {}, Num of edges : {}'.format(G.number_of_nodes(), G.size()))
    emin, emax, edge_weigths = get_minmax_edges(G, 'weight')
    dmin, dmax, edge_length = get_minmax_edges(G, 'length')

    # de = (emax-emin)/100
    # num_edges = []
    # num_nodes = []
    # ener = []
    # for e in xrange(0, 90):
    #     G = filter_edge_by_weight(G, threshold=de*e)
    #     G = remove_node_degree(G) 
    #     num_edges.append(G.size())
    #     num_nodes.append(G.number_of_nodes())
    #     ener.append(de*e)

    # plt.plot(ener, num_nodes)
    # plt.plot(ener, num_edges)
    # plt.show()
    
    # T = nx.minimum_spanning_tree(G)
    # print('Num of nodes : {}, Num of edges : {}'.format(T.number_of_nodes(), T.size()))
    
    # num_bins = 100
    # n, bins, patches = plt.hist(list(nx.get_edge_attributes(G,'weight').values()), num_bins)
    # plt.show()
    
    T = nx.dfs_tree(G, source=1019)
    print(T.edges())
    g = graph()
    
    # C = [G.subgraph(c).copy() for c in sorted(nx.connected_components(G), key=len, reverse=True)]

    # print(' There are {} connected components'.format(len(C)))
    # for i, c in enumerate(C):
    #     if i > 25:
    #         break
    #     print('Num of nodes : {}, Num of edges : {}'.format(c.number_of_nodes(), c.size()))
    #     print_graph(c, "Component_{}.png".format(i))
    
