import matplotlib.pyplot as plt
import numpy as np
import additional as add
import sys
from pathlib import Path
from matplotlib.animation import PillowWriter

fig = plt.figure(figsize=(19, 3.9))
l = plt.plot()

metadata = dict(title='Movie', artist='codinglikemad')
writer = PillowWriter(fps=30, metadata=metadata)

input_path = sys.argv[1]
output_path = sys.argv[2]

#input_path = "D:\\CodeRepos\\Diplom\\Data\\Output\\ToDraw\\2_dim\\Ephi\\"
#output_path = "D:\\CodeRepos\\Diplom\\Drawer\\Pictures\\E_phi\\"

i = 0
path1 = Path(input_path)

with writer.saving(fig, output_path + "EM_field.gif", 100):
    if path1.is_dir():
        for file in path1.iterdir():
            # Read file
            r = []
            z = []
            f = []
            r, z, f = add.read_data(input_path + file.name, r, z, f)
            plt.xlabel('R [m]')
            plt.ylabel('Z [m]')
            plt.title('EM field 2D')

            # Draw mesh
            for ri in r:
                plt.plot([ri, ri], [max(z), min(z)], 'black', linewidth = 0.5)
            for zi in z:
                plt.plot([min(r), max(r)], [zi, zi], 'black', linewidth = 0.5)

            rgrid, zgrid = np.meshgrid(r, z)
            fgrid = []
            i = 0
            it = 0
            for ri in r:
                i = 0
                fgrid.append([])
                for zi in z:
                    fgrid[it].append(float(f[it * len(z) + i]))
                    i+=1
                it+=1

            cbar = plt.contourf(rgrid, zgrid, fgrid, 125, cmap='plasma')
            plt.colorbar(cbar)
            plt.plot([min(r), min(r)], [max(z), min(z)], 'black', linewidth = 3.5)
            plt.plot([min(r), max(r)], [0.0, 0.0], '#b34b06', linewidth = 2.0)
            plt.plot([min(r), max(r)], [-190.0, -190.0], '#351e1c', linewidth = 2.0)
            plt.plot([min(r), max(r)], [-400.0, -400.0], '#efede6', linewidth = 2.0)
            plt.plot([min(r), max(r)], [-812.5, -812.5], '#2c4c3b', linewidth = 2.0)
            name = file.name.replace('.txt', '')
            plt.savefig(output_path + name + '.png')
            i += 1
            #plt.show()
            #input()
            writer.grab_frame()
            plt.clf()
            #plt.close()