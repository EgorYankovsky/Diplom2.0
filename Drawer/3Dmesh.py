# Демонстративная трехмерная сетка.

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


ax = plt.figure(figsize=(19.80, 10.80)).add_subplot(projection='3d')

xarr = [0.0, 1.0, 2.0, 3.0]
yarr = [0.0, 1.0, 2.0, 3.0]
zarr = [0.0, 1.0, 2.0, 3.0]

for k in range(len(zarr)):
    for j in range(len(yarr)):
        for i in range(len(xarr) - 1):
            if xarr[i] == 1.0 and (yarr[j] == 1.0 or yarr[j] == 2.0) and (zarr[k] == 1.0 or zarr[k] == 2.0):
                ax.plot([xarr[i], xarr[i + 1]], [yarr[j], yarr[j]], [zarr[k], zarr[k]], color = "red", linewidth = 3.5)
            elif (j == 0 or j == 3) and (k == 0 or k == 3):
                ax.plot([xarr[i], xarr[i + 1]], [yarr[j], yarr[j]], [zarr[k], zarr[k]], color = "black", linewidth = 2.0)
            else:
                ax.plot([xarr[i], xarr[i + 1]], [yarr[j], yarr[j]], [zarr[k], zarr[k]], color = "black", linewidth = 0.5)

for k in range(len(zarr)):
    for j in range(len(yarr) - 1):
        for i in range(len(xarr)):
            if yarr[j] == 1.0 and (xarr[i] == 1.0 or xarr[i] == 2.0) and (zarr[k] == 1.0 or zarr[k] == 2.0):
                ax.plot([xarr[i], xarr[i]], [yarr[j], yarr[j + 1]], [zarr[k], zarr[k]], color = "red", linewidth = 3.5)
            elif (i == 0 or i == 3) and (k == 0 or k == 3):
                ax.plot([xarr[i], xarr[i]], [yarr[j], yarr[j + 1]], [zarr[k], zarr[k]], color = "black", linewidth = 2.0)
            else:
                ax.plot([xarr[i], xarr[i]], [yarr[j], yarr[j + 1]], [zarr[k], zarr[k]], color = "black", linewidth = 0.5)

for k in range(len(zarr) - 1):
    for j in range(len(yarr)):
        for i in range(len(xarr)):
            if zarr[k] == 1.0 and (yarr[j] == 1.0 or yarr[j] == 2.0) and (xarr[i] == 1.0 or xarr[i] == 2.0):
                ax.plot([xarr[i], xarr[i]], [yarr[j], yarr[j]], [zarr[k], zarr[k + 1]], color = "red", linewidth = 3.5)
            elif (j == 0 or j == 3) and (i == 0 or i == 3):
                ax.plot([xarr[i], xarr[i]], [yarr[j], yarr[j]], [zarr[k], zarr[k + 1]], color = "black", linewidth = 2.0)
            else:
                ax.plot([xarr[i], xarr[i]], [yarr[j], yarr[j]], [zarr[k], zarr[k + 1]], color = "black", linewidth = 0.5)

ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.grid(visible=False)

ax.legend()
ax.set(
    xlabel='X',
    ylabel='Y',
    zlabel='Z')
#plt.savefig('D:\\CodeRepos\\Diplom\\Docs\\Производственная практика\\images\\3D_test_mesh.png')
plt.show()