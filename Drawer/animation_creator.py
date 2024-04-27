import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
import sys
import additional as add
from pathlib import Path

fig = plt.figure(figsize=(19, 10))
l = plt.plot()

plt.xlabel('R [m]')
plt.ylabel('Z [m]')
plt.title('EM field 2D')

metadata = dict(title='Movie', artist='codinglikemad')
writer = PillowWriter(fps=7, metadata=metadata)


relativePath1 = sys.argv[1]
relativePath2 = sys.argv[2]
pictures_path = sys.argv[3]

#relativePath1 = "D:\\CodeRepos\\Diplom\\Data\\Subtotals\\2_dim\\Points.poly"
#relativePath2 = "D:\\CodeRepos\\Diplom\\Data\\Output\\E_phi\\Answer\\"
#pictures_path = "D:\\CodeRepos\\Diplom\\"

add.read_points(relativePath1)

i = 0
path1 = Path(relativePath2)

with writer.saving(fig, pictures_path + "EM_field.gif", 100):
    if path1.is_dir():
        for file in path1.iterdir():
            # Read file
            add.read_basis(relativePath2 + file.name)
            plt.xlabel('R [m]')
            plt.ylabel('Z [m]')
            plt.title('EM field 2D (r,z)')

            # Draw mesh
            for ri in add.r:
                plt.plot([ri, ri], [max(add.z), min(add.z)], 'black', linewidth = 0.5)
            for zi in add.z:
                plt.plot([min(add.r), max(add.r)], [zi, zi], 'black', linewidth = 0.5)
            rgrid, zgrid = np.meshgrid(add.r, add.z)

            #fgrid = f(rgrid, zgrid)
            #cont = plt.contour(rgrid, zgrid, fgrid, 15, colors='black')
            cbar = plt.contourf(rgrid, zgrid, add.basis, 125, cmap='plasma')
            plt.colorbar(cbar)

            plt.plot([min(add.r), min(add.r)], [max(add.z), min(add.z)], 'black', linewidth = 3.5)
            plt.plot([min(add.r), max(add.r)], [0.0, 0.0], 'black', linewidth = 1.8)

            add.basis = []
            writer.grab_frame()
            plt.clf()