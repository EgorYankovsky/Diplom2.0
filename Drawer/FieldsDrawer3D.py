import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection


fig = plt.figure(figsize=(19.80, 10.80))
ax = fig.add_subplot(111, projection='3d')

# воздух
v = np.array([[-5500, -5500, 0],    [ 5500, -5500, 0],    [-5500,  5500, 0],    [ 5500,  5500, 0],
              [-5500, -5500, 1500], [ 5500, -5500, 1500], [-5500,  5500, 1500], [ 5500,  5500, 1500],])
verts = [[v[0], v[1], v[5], v[4]],
         [v[0], v[2], v[6], v[4]]]
ax.add_collection3d(Poly3DCollection(verts, facecolors='#bbeeff', linewidths=1)).set_label("Воздух")

# Глинозём
v = np.array([[-5500, -5500, -500],   [ 5500, -5500, -500],   [-5500,  5500, -500],   [ 5500,  5500, -500],
              [-5500, -5500,  0],     [ 5500, -5500,  0],     [-5500,  5500,  0],     [ 5500,  5500,  0]])
verts = [[v[4], v[5], v[7], v[6]],
         [v[2], v[3], v[7], v[6]],
         [v[1], v[3], v[7], v[5]]]
ax.add_collection3d(Poly3DCollection(verts, facecolors='#b34b06', linewidths=1)).set_label("Слой 1")
# Песчаник
v = np.array([[-5500, -5500, -1000],   [ 5500, -5500, -1000],   [-5500,  5500, -1000],   [ 5500,  5500, -1000],
              [-5500, -5500,  -500],     [ 5500, -5500,  -500],     [-5500,  5500,  -500],     [ 5500,  5500,  -500]])
verts = [[v[4], v[5], v[7], v[6]],
         [v[2], v[3], v[7], v[6]],
         [v[1], v[3], v[7], v[5]]]
ax.add_collection3d(Poly3DCollection(verts, facecolors='#efede6', linewidths=1)).set_label("Слой 2")
# Слюдяной сланец
v = np.array([[-5500, -5500, -1500],   [ 5500, -5500, -1500],   [-5500,  5500, -1500],   [ 5500,  5500, -1500],
              [-5500, -5500,  -1000],     [ 5500, -5500,  -1000],     [-5500,  5500,  -1000],     [ 5500,  5500,  -1000]])
verts = [[v[4], v[5], v[7], v[6]],
         [v[2], v[3], v[7], v[6]],
         [v[1], v[3], v[7], v[5]]]
ax.add_collection3d(Poly3DCollection(verts, facecolors='#2c4c3b', linewidths=1)).set_label("Слой 3")

# серебро
#plt.subplot().add_patch(Rectangle((105, -150.0), 200, 40, facecolor = "#9c9c9c")).set_label("серебро")
# золото
#plt.subplot().add_patch(Rectangle((1205, -180.0), 550, 10, facecolor = "#e5d08f")).set_label("золото")
# грунтовые воды
#plt.subplot().add_patch(Rectangle((2205, -40.0), 150, 10, facecolor = "#3366ff")).set_label("грунтовые воды")
# медь
#plt.subplot().add_patch(Rectangle((955, -85.0), 350, 10, facecolor = "#b27862")).set_label("медь")

ax.plot([-5500.0, -5500.0], [-5500.0,  5500.0], zs = [0.0, 0.0], color = 'black', linewidth=2.0)
ax.plot([-5500.0,  5500.0], [ 5500.0,  5500.0], zs = [0.0, 0.0], color = 'black', linewidth=2.0)
ax.plot([ 5500.0,  5500.0], [ 5500.0, -5500.0], zs = [0.0, 0.0], color = 'black', linewidth=2.0)
ax.plot([ 5500.0, -5500.0], [-5500.0, -5500.0], zs = [0.0, 0.0], color = 'black', linewidth=2.0)

ax.plot([-5500.0, -5500.0], [-5500.0,  5500.0], zs = [1500.0, 1500.0], color = 'black', linewidth=2.0)
ax.plot([-5500.0,  5500.0], [ 5500.0,  5500.0], zs = [1500.0, 1500.0], color = 'black', linewidth=2.0)
ax.plot([ 5500.0,  5500.0], [ 5500.0, -5500.0], zs = [1500.0, 1500.0], color = 'black', linewidth=2.0)
ax.plot([ 5500.0, -5500.0], [-5500.0, -5500.0], zs = [1500.0, 1500.0], color = 'black', linewidth=2.0)

ax.plot([-5500.0, -5500.0], [-5500.0,  5500.0], zs = [-500.0, -500.0], color = 'black', linewidth=2.0)
ax.plot([-5500.0,  5500.0], [ 5500.0,  5500.0], zs = [-500.0, -500.0], color = 'black', linewidth=2.0)
ax.plot([ 5500.0,  5500.0], [ 5500.0, -5500.0], zs = [-500.0, -500.0], color = 'black', linewidth=2.0)
ax.plot([ 5500.0, -5500.0], [-5500.0, -5500.0], zs = [-500.0, -500.0], color = 'black', linewidth=2.0)

ax.plot([-5500.0, -5500.0], [-5500.0,  5500.0], zs = [-1000.0, -1000.0], color = 'black', linewidth=2.0)
ax.plot([-5500.0,  5500.0], [ 5500.0,  5500.0], zs = [-1000.0, -1000.0], color = 'black', linewidth=2.0)
ax.plot([ 5500.0,  5500.0], [ 5500.0, -5500.0], zs = [-1000.0, -1000.0], color = 'black', linewidth=2.0)
ax.plot([ 5500.0, -5500.0], [-5500.0, -5500.0], zs = [-1000.0, -1000.0], color = 'black', linewidth=2.0)

ax.plot([-5500.0, -5500.0], [-5500.0,  5500.0], zs = [-1500.0, -1500.0], color = 'black', linewidth=2.0)
ax.plot([-5500.0,  5500.0], [ 5500.0,  5500.0], zs = [-1500.0, -1500.0], color = 'black', linewidth=2.0)
ax.plot([ 5500.0,  5500.0], [ 5500.0, -5500.0], zs = [-1500.0, -1500.0], color = 'black', linewidth=2.0)
ax.plot([ 5500.0, -5500.0], [-5500.0, -5500.0], zs = [-1500.0, -1500.0], color = 'black', linewidth=2.0)

plt.subplot().legend(loc='upper right')
ax.set(
    xlabel='X [m]',
    ylabel='Y [m]',
    zlabel='Z [m]',
)
plt.show()