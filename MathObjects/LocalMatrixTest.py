import numpy as np
import math

# Config r0, r1, z0, z1
r0 = 1.0
r1 = 1.5
z0 = 1.0
z1 = 1.5

hr = r1 - r0
hz = z1 - z0

Gr = (r0 + 0.5 * hr) / (hr) * np.matrix([[1.0, -1.0],
                                         [-1.0, 1.0]])

Mr = (hr / 6.0) * (r0 * np.matrix([[2.0, 1.0], [1.0, 2.0]]) + 0.5 * hr * np.matrix([[1.0, 1.0], [1.0, 3.0]])) 

d = r0 / hr

Mrr = math.log(1.0 + 1.0 / d) * np.matrix([[(1 + d) ** 2, -1.0 * d * (1.0 + d)], [-1.0 * d * (1.0 + d), d * d]]) - d * np.matrix([[1.0, -1.0], [-1.0, 1.0]]) + 0.5 * np.matrix([[-3.0, 1.0], [1.0, 1.0]])

Gz = (1.0 / hz) * np.matrix([[1.0, -1.0], [-1.0, 1.0]])

Mz = (hz / 6.0) * np.matrix([[2.0, 1.0], [1.0, 2.0]])

for i in range(4):
    for j in range(4):
        print(f"{Gr[i % 2, j % 2] * Mz[i // 2, j // 2] + Mr[i % 2, j % 2] * Gz[i // 2, j // 2] + Mrr[i % 2, j % 2] * Mz[i // 2, j // 2]:.3e}")
    print()


print()
print()

for i in range(4):
    for j in range(4):
        print(f"{Mr[i % 2, j % 2] * Mz[i // 2, j // 2]:.3e}")
    print()