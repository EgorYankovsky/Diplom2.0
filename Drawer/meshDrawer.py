import matplotlib.pyplot as plt
import numpy as np
plt.figure(figsize=(6, 5))

plt.plot([1.0, 2.0], [1.0, 1.0], 'black', linewidth = 1.8)
plt.plot([2.0, 2.0], [1.0, 2.0], 'black', linewidth = 1.8)
plt.plot([1.0, 2.0], [2.0, 2.0], 'black', linewidth = 1.8)
plt.plot([1.0, 1.0], [1.0, 2.0], 'black', linewidth = 1.8)

plt.plot([1.0, 2.0], [4.0 / 3.0, 4.0 / 3.0], 'black', linewidth = 1.1)
plt.plot([1.0, 2.0], [5.0 / 3.0, 5.0 / 3.0], 'black', linewidth = 1.1)

plt.plot([4.0 / 3.0, 4.0 / 3.0], [1.0, 2.0], 'black', linewidth = 1.1)
plt.plot([5.0 / 3.0, 5.0 / 3.0], [1.0, 2.0], 'black', linewidth = 1.1)


plt.show()