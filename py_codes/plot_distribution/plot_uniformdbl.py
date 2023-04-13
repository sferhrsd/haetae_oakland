import matplotlib.pyplot as plt
import random

with open("uniformdbl.txt", "r") as f:
    N = int(f.readline())
    lis = [float(f.readline()) for _ in range(N)]

plt.hist(lis, 100, density=True)


def unif():
    return random.randint(0, 2**32 - 1) * (2 ** (-32))


lis2 = [unif() for _ in range(N)]

plt.show()

plt.hist(lis2, 100, density=True)

plt.show()
