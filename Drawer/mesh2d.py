import matplotlib.pyplot as plt
import numpy as np

r1 = [r for r in range(0, 10)]
r2 = [r for r in range(10, 10111, 111)]

r = r1 + r2

z = [z for z in range(-8000, 2100, 100)]

for ri in r:
    plt.plot([ri, ri], [max(z), min(z)], 'black', linewidth = 0.5)
for zi in z:
    plt.plot([min(r), max(r)], [zi, zi], 'black', linewidth = 0.5)


rgrid, zgrid = np.meshgrid(r, z)


#cont = plt.contour(rgrid, zgrid, fgrid, 15, colors='black')
plt.plot([min(r), max(r)], [0.0, 0.0], '#b34b06', linewidth = 2.0)
plt.plot([0.001, 10000], [-8000.0, -8000.0], color = 'black', linewidth=2.0)
plt.plot([0.001, 0.001], [-8000.0, 2000.0], color = 'black', linewidth=2.0)
plt.plot([10000, 10000], [-8000.0, 2000.0], color = 'black', linewidth=2.0)
plt.plot([0.001, 10000], [2000.0, 2000.0], color = 'grey', linewidth=2.0)
#plt.plot(10.0, 0.0, 'ro') 
#cont.clabel()
plt.legend()
plt.xlabel("R [M]")
plt.ylabel("Z [M]")
plt.show()