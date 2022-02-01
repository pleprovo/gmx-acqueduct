
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter



def potential_6_4(r_on = 1.9, r_off = 4.5, theta_on = 0.0, theta_off = 1.0 , step = 0.01):
    # Make data.
    C = 3855
    D = 738
    R = np.arange(r_on, r_off, step)
    T = np.arange(theta_on, theta_off, step)
    R, T = np.meshgrid(R, R)
    E = (C/np.power(R, 6.0))-(D/np.power(R, 4.0))*np.power(np.cos(T), 4.0)
    return R, T, E


if __name__ == "__main__":
    print "Plop"
    r_on = 1.9
    r_off = 6.0
    theta_on = 0.0
    theta_off = 1.0
    step = 0.05
    R, T, E = potential_6_4(r_on, r_off, theta_on, theta_off, step)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    surf = ax.plot_surface(R, T, E, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    # Customize the z axis.
    # ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    
    plt.show()
