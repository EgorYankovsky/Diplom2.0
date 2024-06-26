import matplotlib.pyplot as plt
import numpy as np
import additional as add
import sys
from pathlib import Path

relativePath1 = sys.argv[1]
relativePath2 = sys.argv[2]
pictures_path = sys.argv[3]

def f(rg, zg): return ((rg - 10) ** 2 + zg ** 2) ** 0.5

add.read_points(relativePath1)

i = 0
path1 = Path(relativePath2)

if path1.is_dir():
    for file in path1.iterdir():
        if file.name != "Discrepancy.dat":
            # Read file
            add.read_basis(relativePath2 + file.name)

            # Generate plot
            plt.figure(figsize=(19, 10))

            # Draw (r,z) axis
            #plt.plot([min(add.r), min(add.r)], [max(add.z), min(add.z)], 'black', linewidth = 3.5)
            #plt.plot([min(add.r), max(add.r)], [0.0, 0.0], 'black', linewidth = 1.8)
            plt.plot([min(add.r), max(add.r)], [-60.0, -60.0], 'black', linewidth = 1.0)
            plt.plot([min(add.r), max(add.r)], [-40.0, -40.0], 'black', linewidth = 1.0)
            plt.plot([min(add.r), max(add.r)], [-20.0, -20.0], 'black', linewidth = 1.0)

            # Draw mesh
            for ri in add.r:
                plt.plot([ri, ri], [max(add.z), min(add.z)], 'black', linewidth = 0.5)
            for zi in add.z:
                plt.plot([min(add.r), max(add.r)], [zi, zi], 'black', linewidth = 0.5)


            rgrid, zgrid = np.meshgrid(add.r, add.z)
            fgrid = f(rgrid, zgrid)

            cont = plt.contour(rgrid, zgrid, fgrid, 15, colors='black')
            cbar = plt.contourf(rgrid, zgrid, add.basis, 125, cmap='plasma')
            plt.colorbar(cbar)
            plt.title(file.name)
            plt.plot([min(add.r), min(add.r)], [max(add.z), min(add.z)], 'black', linewidth = 3.5)
            plt.plot([min(add.r), max(add.r)], [0.0, 0.0], 'black', linewidth = 1.8)
            plt.plot(10.0, 0.0, 'ro') 
            cont.clabel()

            name = file.name.replace('.dat', '')

            plt.savefig(pictures_path + name + '.png')
            i += 1
            add.basis = []
            #plt.show()
            plt.close()