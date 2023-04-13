import matplotlib.pyplot as plt
import numpy as np

dim = 256 * (3 + 3)  # dimension of hyperball

B = 10705.083778907292  # B(radius of hyperball)

with open("hyperball_radius.txt", "r") as f:
    N = int(f.readline())
    lis = [float(f.readline()) for _ in range(N)]

plt.hist(lis, 100, density=True)

xplt = np.linspace(0, B, 10000)
yplt = dim * ((xplt / B) ** (dim - 1)) / B
# PDF : n \cdot x^{n - 1}
plt.plot(xplt, yplt)
plt.xlim(0, B)
plt.show()

plt.hist(lis, 100, density=True)
plt.plot(xplt, yplt)
plt.xlim(B - B // 100, B)  # 99% ~ 100%
plt.show()

plt.hist(lis, 100, density=True)
plt.plot(xplt, yplt)
plt.xlim(B - B // 1000, B)  # 99.9% ~ 100%
plt.show()
