import random
chset = "0123456789abcdef"
def f():
    return "0x" + "".join([random.choice(chset) for _ in range(8)])
lis = [f() for _ in range(24)]
print("{", end="")
print(*lis, end="", sep=", ")
print("}")
