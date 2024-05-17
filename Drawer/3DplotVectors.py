# Отрисовывает вектор-функцию

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import additional as add
from matplotlib.animation import PillowWriter
from pathlib import Path
import sys

#read_path = sys.argv[1]
#write_path = sys.argv[2]

read_path = "D:\\CodeRepos\\Diplom\\Data\\Output\\E_phi\\ToDraw\\ConvertedTo3D\\"
write_path = "D:\\CodeRepos\\Diplom\\Drawer\\Pictures\\E_phi3D\\"

#metadata = dict(title='Movie', artist='codinglikemad')
#writer = PillowWriter(fps=7, metadata=metadata)

read_dir = Path(read_path)
write_dir = Path(write_path)

ax = plt.figure(figsize=(19.80, 10.80)).add_subplot(projection='3d')

#with writer.saving(ax, write_path + "EM_field.gif", 100):
if read_dir.is_dir():
    for file in read_dir.iterdir():
        add.read_vectors(read_path + file.name)
        ax.set_xlim(add.max_values[0], add.max_values[1])
        ax.set_ylim(add.max_values[2], add.max_values[3])
        ax.set_zlim(add.max_values[4], add.max_values[5])
        for rib in add.ribs_arr:
            coord = rib.get_values_to_draw()
            ax.quiver(coord[0], coord[1], coord[2], coord[3], coord[4], coord[5], length = 50.0, normalize=True, color='red')
        ax.set(
            xlabel='X',
            ylabel='Y',
            zlabel='Z')
        name = file.name.replace('.txt', '')
        plt.savefig(write_path + name + '.png')
        #plt.show()