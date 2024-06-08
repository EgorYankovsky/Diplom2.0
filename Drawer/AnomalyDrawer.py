import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

plt.figure(figsize=(9.9, 3.9))
plt.plot([-5500.0, 5500.0], [-9000.0, -9000.0], color = 'black', linewidth=1.0)
plt.plot([-5500.0, -5500.0], [-9000.0, 3000.0], color = 'black', linewidth=1.0)
plt.plot([5500.0, 5500.0], [-9000.0, 3000.0], color = 'black', linewidth=1.0)
plt.plot([-5500.0, 5500.0], [3000.0, 3000.0], color = 'black', linewidth=1.0)

plt.plot([-5500.0, 5500.0], [0.0, 0.0], color = 'black', linewidth=1.0)
plt.plot([-5500.0, 5500.0], [-3000.0, -3000.0], color = 'black', linewidth=1.0)
plt.plot([-5500.0, 5500.0], [-6000.0, -6000.0], color = 'black', linewidth=1.0)

# воздух
plt.subplot().add_patch(Rectangle((-5500.0, 0.0), 11000.0, 3000, facecolor = "#bbeeff")).set_label("Воздух")
# Глинозём
plt.subplot().add_patch(Rectangle((-5500.0, -3000.0), 11000.0, 3000, facecolor = "#b34b06")).set_label("Слой 1")
# Песчаник
plt.subplot().add_patch(Rectangle((-5500.0, -6000.0), 11000.0, 3000, facecolor = "#efede6")).set_label("Слой 2")
# Слюдяной сланец
plt.subplot().add_patch(Rectangle((-5500.0, -9000.0), 11000.0, 3000, facecolor = "#2c4c3b")).set_label("Слой 3")
# серебро
plt.subplot().add_patch(Rectangle((1000.0, -2500.0), 1500.0, 1500, facecolor = "#9c9c9c")).set_label("Аномалия 1")
# золото
plt.subplot().add_patch(Rectangle((-3500.0, -3500.0), 1500.0, 1500, facecolor = "#e5d08f")).set_label("Аномалия 2")


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
    xlabel='Y [M]',
    ylabel='Z [M]',)
plt.show()