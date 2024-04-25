import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

plt.figure(figsize=(19.80, 10.80))
fig, ax = plt.subplots()
#ax.plot([0.001, 10000], [0.0, 0.0], color = 'black', linewidth=2.0)

ax.plot([0.001, 10000], [-8000.0, -8000.0], color = 'black', linewidth=2.0)
ax.plot([0.001, 0.001], [-8000.0, 2000.0], color = 'black', linewidth=2.0)
ax.plot([10000, 10000], [-8000.0, 2000.0], color = 'black', linewidth=2.0)
ax.plot([0.001, 10000], [2000.0, 2000.0], color = 'grey', linewidth=2.0)
ax.plot(10, 0.0, 'ro', color = 'red')


ax.add_patch(Rectangle((0.001, -8000.0), 9999.999, 8000, facecolor = "#b34b06"))
ax.add_patch(Rectangle((0.001, 0.0), 9999.999, 2000, facecolor = "#1f4af7"))

ax.legend()
ax.set(
    xlabel='R [M]',
    ylabel='Z [M]',)

plt.show()