
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm

fig = plt.figure()
ax = fig.add_subplot(111, projection='polar')

# Make data.
C = 3855
D = 738
ron = 5.5
roff = 6.5
ton = 0.25
toff = 0.0301

# Create the mesh in polar coordinates and compute corresponding Z.
r = np.linspace(2.5, 7.0, 100)
p = np.linspace(-np.pi/2, np.pi/2, 100)
pb = np.cos(p)
P, R = np.meshgrid(pb, r)
Z = ((C/pow(R, 6.0))-(D/pow(R, 4.0)))*pow(P, 4.0)
print Z.shape
# Express the mesh in the cartesian system.
# X, Y = R*np.cos(P), R*np.sin(P)

# Plot the surface.
ctf = ax.contourf(p, r, Z, 100, cmap=cm.jet)
plt.colorbar(ctf)

# Tweak the limits and add latex math labels.
#ax.set_zlim(-5, 5)
#ax.set_xlabel(r'$\theta_\mathrm{real}$')
#ax.set_ylabel(r'$\theta_\mathrm{im}$')
#ax.set_zlabel(r'$V(\theta)$')

plt.show()



