import math

k = 2
l = 4
n = 1 << 8
alpha = 1 << 8


def computeB1(B2):
    return B2 - math.sqrt(n * (k + l)) / 2 - alpha * math.sqrt(n * k) / 4


def computeB0(B1):
    exponent = 1 / (n * (k + l))
    term = math.pow(2, exponent)
    return B1 * term


q = [67073, 57857, 65025, 45569, 88577]
B2 = [12607.212, 12607.212, 11190.369619213772, 16144.820319078379, 23693.57972306502]
for i in range(0, 5):
    B1 = computeB1(B2[i])
    B0 = computeB0(B1)
    print(
        "For q = {0}, (B, B', B'') = ({1:18.12f}, {2:18.12f}, {3:18.12f})".format(
            q[i], B0, B1, B2[i]
        )
    )
