
from numpy import asarray
from mayavi.mlab import points3d, quiver3d

def read_triangulation(file):
    lines = [ line.split() for line in file.readlines() ]
    points = asarray(lines[2:int(lines[0][0])+2]).astype(float)
    plop = asarray(lines[int(lines[0][0])+2:-2])
    edges = plop[:,:-1].astype(int)
    weights = plop[:,2].astype(float)
    source = asarray(lines[-2]).astype(float)
    sink = asarray(lines[-1]).astype(float)
    return points, edges, weights, source, sink

def plot_mayavi(points, edges, weights, source, sink):
    xq = []
    yq = []
    zq = []
    uq = []
    vq = []
    wq = []
    weight_max = max(weights)
    for i, edge in enumerate(edges):
        xq.append(points[edge[0]][0])
        yq.append(points[edge[0]][1])
        zq.append(points[edge[0]][2])
        uq.append(points[edge[1]][0]-points[edge[0]][0])
        vq.append(points[edge[1]][1]-points[edge[0]][1])
        wq.append(points[edge[1]][2]-points[edge[0]][2])
    point = points3d(xq, yq, zq, scale_factor=0.1, color=(1.0, 0.0, 0.0))
    #sources = points3d(source[0], source[1], source[2], scale_factor=0.1, color=(0.0, 1.0, 0.0))
    #sinks = points3d(sink[0], sink[1], sink[2], scale_factor=0.1, color=(0.0, 0.0, 1.0))
    flow = quiver3d(xq, yq, zq, uq, vq, wq, scale_mode='scalar', scalars=weights, vmax=weight_max, vmin=0.0, mode='arrow')
    #mlab.outline()
    #mlab.show()
    return point, flow, sources, sinks

