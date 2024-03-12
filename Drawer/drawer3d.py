import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(19.80, 10.80))

ax_3d = fig.add_subplot(projection='3d')

ax_3d.set_xlabel('r')
ax_3d.set_ylabel('z')
ax_3d.set_zlabel('u')

r = np.arange(0, 10, 1.0)
z = np.arange(0, 10, 1.0)

rgrid, zgrid = np.meshgrid(r, z)

ugrid = rgrid ** 2 * zgrid ** 2
ugrid1 = 1.1 * rgrid ** 2 * zgrid ** 2

ax_3d.plot_surface(rgrid, zgrid, ugrid, cmap='plasma')
ax_3d.plot_wireframe(rgrid, zgrid, ugrid1)

plt.show()