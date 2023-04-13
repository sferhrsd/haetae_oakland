import matplotlib.pyplot as plt
import numpy as np
import math

with open("gaussian.txt", "r") as f:
    N = int(f.readline())
    lis = [float(f.readline()) for _ in range(N)]

plt.hist(lis, 100, density=True)

xplt = np.linspace(-3, 3, 100)
yplt = 1 / math.sqrt(2 *  math.pi) * np.exp(-xplt * xplt / 2)
plt.plot(xplt, yplt)
plt.show()
