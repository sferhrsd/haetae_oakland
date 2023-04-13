import math

N = 256


def f(x):
    return math.exp(-x * x / 2)


def invf(y):
    return math.sqrt(-2 * math.log(y))


s = 1.0
e = 6.0
for rep in range(100000):
    flag = False
    m = (s + e) / 2
    x = [0, m]
    y = [0]
    t = math.sqrt(math.pi / 2) * math.erfc(x[1] / math.sqrt(2))
    y.append(f(x[1]))
    A = x[1] * y[1] + t
    x[0] = A / y[1]
    while len(x) < N + 1:
        y.append(y[-1] + A / x[-1])
        if y[-1] > f(0):
            flag = True
            break
        x.append(invf(y[-1]))
    if flag:
        s = m
    else:
        e = m
print(s)
print(e)
