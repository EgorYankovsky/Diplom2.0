# Отрисовывает вектор-функцию

import matplotlib.pyplot as plt
import numpy as np

ax = plt.figure(figsize=(19.80, 10.80)).add_subplot(projection='3d')

# Make the grid
x, y, z = np.meshgrid(np.arange(start=-8, stop=10, step=3),
                      np.arange(start=-8, stop=10, step=3),
                      np.arange(start=-8, stop=10, step=3))

# Make the direction data for the arrows
Ax = np.sin(z) * np.cos(y)
Ay = np.sin(x) * np.cos(z)
Az = np.sin(x) * np.cos(y)

#Ex = z
#Ey = x
#Ez = y

ax.quiver(x, y, z, Ax, Ay, Az, length=1, color='blue')
#ax.quiver(x, y, z, Ex, Ey, Ez, length=0.2, color='red')
ax.legend()
ax.set(
    xlabel='X',
    ylabel='Y',
    zlabel='Z')
plt.show()