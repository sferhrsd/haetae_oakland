import math

N = 256


def f(x):
    return math.exp(-x * x / 2)


def invf(y):
    return math.sqrt(-2 * math.log(y))


x1 = 3.6541528853610091
x = [0, x1]
y = [0]
t = math.sqrt(math.pi / 2) * math.erfc(x[1] / math.sqrt(2))
y.append(f(x[1]))
A = x[1] * y[1] + t
x[0] = A / y[1]
while len(x) < N + 1:
    y.append(y[-1] + A / x[-1])
    y[-1] = min(y[-1], 1.0)
    x.append(invf(y[-1]))
k = []
for i in range(len(x)):
    if i != len(x) - 1:
        k.append(math.ceil(x[i + 1] / x[i] * (2**31)))
    x[i] *= 2 ** (-31)

x = x[:-1]
y = y[1:]

print(x)
print(y)
print(k)


def print_c(lis):
    print("{", end="")
    print(*lis, end="", sep=", ")
    print("}")


#print_c(x)
#print_c(y)
#print_c(k)
