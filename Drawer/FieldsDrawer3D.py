import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection


fig = plt.figure(figsize=(19.80, 10.80))
ax = fig.add_subplot(111, projection='3d')

# воздух
v = np.array([[-2500, -2500, 0],   [ 2500, -2500, 0],   [-2500,  2500, 0],   [ 2500,  2500, 0],
              [-2500, -2500, 200], [ 2500, -2500, 200], [-2500,  2500, 200], [ 2500,  2500, 200],])
verts = [[v[0], v[1], v[5], v[4]],
         [v[0], v[2], v[6], v[4]]]
ax.add_collection3d(Poly3DCollection(verts, facecolors='#bbeeff', linewidths=1)).set_label("Воздух")

# Глинозём
v = np.array([[-2500, -2500, -9],   [ 2500, -2500, -9],   [-2500,  2500, -9],   [ 2500,  2500, -9],
              [-2500, -2500,  0],   [ 2500, -2500,  0],   [-2500,  2500,  0],   [ 2500,  2500,  0]])
verts = [[v[4], v[5], v[7], v[6]],
         [v[2], v[3], v[7], v[6]],
         [v[1], v[3], v[7], v[5]]]
ax.add_collection3d(Poly3DCollection(verts, facecolors='#b34b06', linewidths=1)).set_label("Глинозём")
# Песчаник
v = np.array([[-2500, -2500, -50],   [ 2500, -2500, -50],   [-2500,  2500, -50],   [ 2500,  2500, -50],
              [-2500, -2500,  -9],   [ 2500, -2500,  -9],   [-2500,  2500,  -9],   [ 2500,  2500,  -9]])
verts = [[v[4], v[5], v[7], v[6]],
         [v[2], v[3], v[7], v[6]],
         [v[1], v[3], v[7], v[5]]]
ax.add_collection3d(Poly3DCollection(verts, facecolors='#efede6', linewidths=1)).set_label("Песчаник")
# Известняк, юрский
v = np.array([[-2500, -2500, -125],   [ 2500, -2500, -125],   [-2500,  2500, -125],   [ 2500,  2500, -125],
              [-2500, -2500,  -50],   [ 2500, -2500,  -50],   [-2500,  2500,  -50],   [ 2500,  2500,  -50]])
verts = [[v[4], v[5], v[7], v[6]],
         [v[2], v[3], v[7], v[6]],
         [v[1], v[3], v[7], v[5]]]
ax.add_collection3d(Poly3DCollection(verts, facecolors='#351e1c', linewidths=1)).set_label("Известняк, юрский")
# Слюдяной сланец
v = np.array([[-2500, -2500, -200],   [ 2500, -2500, -200],   [-2500,  2500, -200],   [ 2500,  2500, -200],
              [-2500, -2500, -125],   [ 2500, -2500, -125],   [-2500,  2500, -125],   [ 2500,  2500, -125]])
verts = [[v[4], v[5], v[7], v[6]],
         [v[2], v[3], v[7], v[6]],
         [v[1], v[3], v[7], v[5]]]
ax.add_collection3d(Poly3DCollection(verts, facecolors='#2c4c3b', linewidths=1)).set_label("Слюдяной сланец")

# серебро
#plt.subplot().add_patch(Rectangle((105, -150.0), 200, 40, facecolor = "#9c9c9c")).set_label("серебро")
# золото
#plt.subplot().add_patch(Rectangle((1205, -180.0), 550, 10, facecolor = "#e5d08f")).set_label("золото")
# грунтовые воды
#plt.subplot().add_patch(Rectangle((2205, -40.0), 150, 10, facecolor = "#3366ff")).set_label("грунтовые воды")
# медь
#plt.subplot().add_patch(Rectangle((955, -85.0), 350, 10, facecolor = "#b27862")).set_label("медь")

ax.plot([-2500.0, -2500.0], [-2500.0,  2500.0], zs = [0.0, 0.0], color = 'black', linewidth=2.0)
ax.plot([-2500.0,  2500.0], [ 2500.0,  2500.0], zs = [0.0, 0.0], color = 'black', linewidth=2.0)
ax.plot([ 2500.0,  2500.0], [ 2500.0, -2500.0], zs = [0.0, 0.0], color = 'black', linewidth=2.0)
ax.plot([ 2500.0, -2500.0], [-2500.0, -2500.0], zs = [0.0, 0.0], color = 'black', linewidth=2.0)

ax.plot([-2500.0, -2500.0], [-2500.0,  2500.0], zs = [200.0, 200.0], color = 'black', linewidth=2.0)
ax.plot([-2500.0,  2500.0], [ 2500.0,  2500.0], zs = [200.0, 200.0], color = 'black', linewidth=2.0)
ax.plot([ 2500.0,  2500.0], [ 2500.0, -2500.0], zs = [200.0, 200.0], color = 'black', linewidth=2.0)
ax.plot([ 2500.0, -2500.0], [-2500.0, -2500.0], zs = [200.0, 200.0], color = 'black', linewidth=2.0)

ax.plot([-2500.0, -2500.0], [-2500.0,  2500.0], zs = [-9.0, -9.0], color = 'black', linewidth=2.0)
ax.plot([-2500.0,  2500.0], [ 2500.0,  2500.0], zs = [-9.0, -9.0], color = 'black', linewidth=2.0)
ax.plot([ 2500.0,  2500.0], [ 2500.0, -2500.0], zs = [-9.0, -9.0], color = 'black', linewidth=2.0)
ax.plot([ 2500.0, -2500.0], [-2500.0, -2500.0], zs = [-9.0, -9.0], color = 'black', linewidth=2.0)

ax.plot([-2500.0, -2500.0], [-2500.0,  2500.0], zs = [-50.0, -50.0], color = 'black', linewidth=2.0)
ax.plot([-2500.0,  2500.0], [ 2500.0,  2500.0], zs = [-50.0, -50.0], color = 'black', linewidth=2.0)
ax.plot([ 2500.0,  2500.0], [ 2500.0, -2500.0], zs = [-50.0, -50.0], color = 'black', linewidth=2.0)
ax.plot([ 2500.0, -2500.0], [-2500.0, -2500.0], zs = [-50.0, -50.0], color = 'black', linewidth=2.0)

ax.plot([-2500.0, -2500.0], [-2500.0,  2500.0], zs = [-125.0, -125.0], color = 'black', linewidth=2.0)
ax.plot([-2500.0,  2500.0], [ 2500.0,  2500.0], zs = [-125.0, -125.0], color = 'black', linewidth=2.0)
ax.plot([ 2500.0,  2500.0], [ 2500.0, -2500.0], zs = [-125.0, -125.0], color = 'black', linewidth=2.0)
ax.plot([ 2500.0, -2500.0], [-2500.0, -2500.0], zs = [-125.0, -125.0], color = 'black', linewidth=2.0)

ax.plot([-2500.0, -2500.0], [-2500.0,  2500.0], zs = [-200.0, -200.0], color = 'black', linewidth=2.0)
ax.plot([-2500.0,  2500.0], [ 2500.0,  2500.0], zs = [-200.0, -200.0], color = 'black', linewidth=2.0)
ax.plot([ 2500.0,  2500.0], [ 2500.0, -2500.0], zs = [-200.0, -200.0], color = 'black', linewidth=2.0)
ax.plot([ 2500.0, -2500.0], [-2500.0, -2500.0], zs = [-200.0, -200.0], color = 'black', linewidth=2.0)

ax.plot([-2500.0, -2500.0], [-2500.0, -2500.0], zs = [-200.0, 200.0], color = 'black', linewidth=2.0)
ax.plot([-2500.0, -2500.0], [ 2500.0,  2500.0], zs = [-200.0, 200.0], color = 'black', linewidth=2.0)
ax.plot([ 2500.0,  2500.0], [-2500.0, -2500.0], zs = [-200.0, 200.0], color = 'black', linewidth=2.0)
ax.plot([ 2500.0,  2500.0], [ 2500.0,  2500.0], zs = [-200.0, 200.0], color = 'black', linewidth=2.0)

plt.subplot().legend(loc='upper right')
ax.set(
    xlabel='X [m]',
    ylabel='Y [m]',
    zlabel='Z [m]',
)
plt.show()