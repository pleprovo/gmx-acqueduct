from network_vis import *
from alphashape_vis import *
from mayavi.mlab import show


if __name__ == "__main__":
    with open('tri.dat') as file:
        vertices, faces = read_off(file)

    plot_mayavi_off(vertices, faces)

    show()
