import matplotlib.pyplot as plt
import numpy as np

r0 = 10.0
hr = 20.14
r = [r0 + hr * i for i in range(176)]
r.insert(0, 0)
z0 = -200.0
hz = 6.67
z = [z0 + hz * i for i in range(61)]

for ri in r:
    plt.plot([ri, ri], [max(z), min(z)], 'black', linewidth = 0.5)
for zi in z:
    plt.plot([min(r), max(r)], [zi, zi], 'black', linewidth = 0.5)

rgrid, zgrid = np.meshgrid(r, z)


#cont = plt.contour(rgrid, zgrid, fgrid, 15, colors='black')
plt.plot([min(r), max(r)], [0.0, 0.0], '#b34b06', linewidth = 2.0)
plt.plot([min(r), max(r)], [-9.0, -9.0], '#efede6', linewidth = 2.0)
plt.plot([min(r), max(r)], [-50.0, -50.0], '#351e1c', linewidth = 2.0)
plt.plot([min(r), max(r)], [-125.0, -125.0], '#2c4c3b', linewidth = 2.0)

plt.plot([0.001, 3535], [-200.0, -200.0], color = 'black', linewidth=2.0)
plt.plot([0.001, 0.001], [-200.0, 200.0], color = 'black', linewidth=2.0)
plt.plot([3535, 3535], [-200.0, 200.0], color = 'black', linewidth=2.0)
plt.plot([0.001, 3535], [200.0, 200.0], color = 'black', linewidth=2.0)
#plt.plot(10.0, 0.0, 'ro') 
#cont.clabel()
plt.legend()
plt.xlabel("R [M]")
plt.ylabel("Z [M]")
plt.show()