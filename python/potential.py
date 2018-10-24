
from numpy import asarray, loadtxt
import matplotlib.pyplot as plt

plt.figure(1)
mat = loadtxt("potential.dat", dtype=float)   
plt.imshow(mat)

plt.figure(2)  
radius = loadtxt("switch_radius.dat",dtype=float)
plt.plot(radius[0,:], radius[1,:])

angle = loadtxt("switch_angle.dat",dtype=float)
plt.plot(angle[0,:], angle[1,:])
plt.show()
