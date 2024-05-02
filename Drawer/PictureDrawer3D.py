import matplotlib.pyplot as plt
import numpy as np
import additional as add
import sys
from pathlib import Path

#relativePath1 = sys.argv[1]
#relativePath2 = sys.argv[2]
#pictures_path = sys.argv[3]

points_path = "D:\\CodeRepos\\Diplom\\Data\\Subtotals\\3_dim\\Field0\\Points.poly"
ribs_path = "D:\\CodeRepos\\Diplom\\Data\\Subtotals\\3_dim\\Field0\\Ribs.poly"
answer_path = "D:\\CodeRepos\\Diplom\\Data\\Output\\A_phi\\Answer3D\\ConvertedTo3D\\"
pictures_path = "D:\\CodeRepos\\Diplom\\"

add.read_ribs(ribs_path)
add.read_points(points_path)

i = 0
path1 = Path(answer_path)

if path1.is_dir():
    for file in path1.iterdir():
        add.read_basis(answer_path + file.name)

        plt.figure(figsize=(19, 10), projection = '3d')

        