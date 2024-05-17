import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

plt.figure(figsize=(9.9, 2.9))
plt.plot([0.001, 3535.533], [-200.0, -200.0], color = 'black', linewidth=1.0)
plt.plot([0.001, 0.001], [-200.0, 200.0], color = 'black', linewidth=1.0)
plt.plot([3535.533, 3535.533], [-200.0, 200.0], color = 'black', linewidth=1.0)
plt.plot([0.001, 3535.533], [200.0, 200.0], color = 'black', linewidth=1.0)

plt.plot([0.001, 3535.533], [0.0, 0.0], color = 'black', linewidth=1.0)
plt.plot([0.001, 3535.533], [-9.0, -9.0], color = 'black', linewidth=1.0)
plt.plot([0.001, 3535.533], [-50.0, -50.0], color = 'black', linewidth=1.0)
plt.plot([0.001, 3535.533], [-125.0, -125.0], color = 'black', linewidth=1.0)


# воздух
plt.subplot().add_patch(Rectangle((0.001, 0.0), 3535.532, 200, facecolor = "#bbeeff")).set_label("Воздух")
# Глинозём
plt.subplot().add_patch(Rectangle((0.001, -9.0), 3535.532, 9, facecolor = "#b34b06")).set_label("Глинозём")
# Песчаник
plt.subplot().add_patch(Rectangle((0.001, -50.0), 3535.532, 41, facecolor = "#efede6")).set_label("Песчаник")
# Известняк, юрский
plt.subplot().add_patch(Rectangle((0.001, -125.0), 3535.532, 75, facecolor = "#351e1c")).set_label("Юрский известняк")
# Слюдяной сланец
plt.subplot().add_patch(Rectangle((0.001, -200.0), 3535.532, 75, facecolor = "#2c4c3b")).set_label("Слюдяной сланец")

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