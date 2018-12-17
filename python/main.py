
from network_vis import *
from alphashape_vis import *
from mayavi.mlab import show
#import matplotlib as mpl
#from mpl_toolkits.mplot3d import Axes3D

if __name__ == "__main__":
    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    with open('triangulation.dat') as file:
        points, edges, weights, source, sink = read_triangulation(file)
    
    with open('surface.off') as file:
         vertices, faces = read_off(file)
         
    plot_mayavi_off(vertices, faces)
    plot_mayavi(points, edges, weights, source, sink)
    
    show()
