import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

plt.figure(figsize=(9.9, 2.9))
plt.plot([0.001, 77781.75], [-9000.0, -9000.0], color = 'black', linewidth=1.0)
plt.plot([0.001, 0.001], [-9000.0, 3000.0], color = 'black', linewidth=1.0)
plt.plot([77781.75, 77781.75], [-9000.0, 3000.0], color = 'black', linewidth=1.0)
plt.plot([0.001, 77781.75], [3000.0, 3000.0], color = 'black', linewidth=1.0)

plt.plot([0.001, 77781.75], [0.0, 0.0], color = 'black', linewidth=1.0)
plt.plot([0.001, 77781.75], [-3000.0, -3000.0], color = 'black', linewidth=1.0)
plt.plot([0.001, 77781.75], [-6000.0, -6000.0], color = 'black', linewidth=1.0)


# воздух
plt.subplot().add_patch(Rectangle((0.001, 0.0), 77781.75, 3000, facecolor = "#bbeeff")).set_label("Воздух")
# Глинозём
plt.subplot().add_patch(Rectangle((0.001, -3000.0), 77781.75, 3000, facecolor = "#b34b06")).set_label("Слой 1")
# Песчаник
plt.subplot().add_patch(Rectangle((0.001, -6000.0), 77781.75, 3000, facecolor = "#efede6")).set_label("Слой 2")
# Слюдяной сланец
plt.subplot().add_patch(Rectangle((0.001, -9000.0), 77781.75, 3000, facecolor = "#2c4c3b")).set_label("Слой 3")

# серебро
#plt.subplot().add_patch(Rectangle((105, -150.0), 200, 40, facecolor = "#9c9c9c")).set_label("серебро")
# золото
#plt.subplot().add_patch(Rectangle((1205, -180.0), 550, 10, facecolor = "#e5d08f")).set_label("золото")
# грунтовые воды
#plt.subplot().add_patch(Rectangle((2205, -40.0), 150, 10, facecolor = "#3366ff")).set_label("грунтовые воды")
# медь
#plt.subplot().add_patch(Rectangle((955, -85.0), 350, 10, facecolor = "#b27862")).set_label("медь")
plt.plot()
plt.subplot().legend(loc='upper right')
plt.subplot().set(
    xlabel='R [M]',
    ylabel='Z [M]',)
plt.show()