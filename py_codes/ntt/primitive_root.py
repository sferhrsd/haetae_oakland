q = 12289 # prime for ring
N = 256
# finds x^2N \equiv 1 \pmod q , x^i \not\equiv 1 \pmod q (1 \leq i < 2N)

lis = []

for x in range(q):
    if (x ** (2 * N)) % q != 1:
        continue
    flag = True
    now = x
    for i in range(1, 2 * N):
        if now == 1:
           flag = False
           break
        now = (now * x) % q
    if flag:
        lis.append(x)

print(len(lis))
print(lis)

