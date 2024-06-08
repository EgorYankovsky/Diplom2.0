import matplotlib.pyplot as plt
import numpy as np
import additional as add
import sys
import os
from pathlib import Path
from matplotlib.animation import PillowWriter

fig = plt.figure(figsize=(19, 3.9))
l = plt.plot()

metadata = dict(title='Movie', artist='codinglikemad')
writer = PillowWriter(fps=12, metadata=metadata)

#input_path = sys.argv[1]
#output_path = sys.argv[2]

#input_path = "D:\\CodeRepos\\Diplom\\Data\\Output\\ToDraw\\2_dim\\Ephi\\"
#output_path = "D:\\CodeRepos\\Diplom\\Drawer\\Pictures\\E_phi\\"

input_path = "D:\\CodeRepos\\Diplom\\Data\\Output\\ToDraw\\3_dim\\E\\"
output_path = "D:\\CodeRepos\\Diplom\\Drawer\\Pictures\\E3d\\"

files_amount = len(os.listdir(input_path))

i = 0
current_file_num = 0
path1 = Path(input_path)
with writer.saving(fig, output_path + "EM_field.gif", 100):
    if path1.is_dir():
        for current_file_num in range(files_amount):
            for file in path1.iterdir():
                # Read file
                if (file.name != f'Answer_E_time_layer_{current_file_num}.txt'):
                    continue
                r = []
                z = []
                f = []
                r, z, f = add.read_data(input_path + file.name, r, z, f)
                plt.xlabel('R [m]')
                plt.ylabel('Z [m]')
                plt.title('EM field 2D')

                # Draw mesh
                #for ri in r:
                #    plt.plot([ri, ri], [max(z), min(z)], 'black', linewidth = 0.5)
                #for zi in z:
                #    plt.plot([min(r), max(r)], [zi, zi], 'black', linewidth = 0.5)

                #plt.plot([2548, 2548], [-180, -80], 'black', linewidth=2.0)
                #plt.plot([2731, 2731], [-180, -80], 'black', linewidth=2.0)
                #plt.plot([2548, 2731], [-80, -80], 'black', linewidth=2.0)
                #plt.plot([2548, 2731], [-180, -180], 'black', linewidth=2.0)

                #plt.plot([-2605, -2605], [-1250, -800], 'black', linewidth=2.0)
                #plt.plot([-2415, -2415], [-1250, -800], 'black', linewidth=2.0)
                #plt.plot([-2605, -2415], [-800, -800], 'black', linewidth=2.0)
                #plt.plot([-2605, -2415], [-1250, -1250], 'black', linewidth=2.0)

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


                cbar = plt.contourf(rgrid, zgrid, fgrid, 125, cmap='jet')
                plt.colorbar(cbar)
                plt.plot([min(r), min(r)], [max(z), min(z)], 'black', linewidth = 3.5)
                plt.plot([min(r), max(r)], [0.0, 0.0], '#b34b06', linewidth = 2.0)
                plt.plot([min(r), max(r)], [-3000.0, -3000.0], '#351e1c', linewidth = 2.0)
                plt.plot([min(r), max(r)], [-6000.0, -6000.0], '#efede6', linewidth = 2.0)
                plt.scatter(1000, 0, color='green')
                plt.scatter(100, 0, color='blue')
                plt.scatter(2000, 0, color='pink')
                plt.scatter(5000, 0, color='red')
                name = file.name.replace('.txt', '')
                plt.savefig(output_path + name + '.png')
                i += 1
                #plt.show()
                #input()
                writer.grab_frame()
                plt.clf()
                #plt.close()