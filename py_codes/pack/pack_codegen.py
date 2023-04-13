from math import gcd

s = -10700  # int(input())
e = 10700  # int(input())
ename = "10701"  # input()

blen = (e - s).bit_length()
print(f"bit length : {blen}")
typeint = [32, 16, 8]
typestr = ""
for num in typeint:
    if num >= blen:
        typestr = f"uint{num}_t"
print("typestr: ", typestr)
andstr = hex(2**blen - 1)
print("andstr: ", andstr)


g = gcd(blen, 8)
tcnt = 8 // g
rcnt = blen // g
print(f"gcd : {g} , itercnt : {tcnt}")


print("pack part")


def getj(j, mode):
    if mode:
        return f"t[{j}]"
    else:
        return f"a->coeffs[{tcnt}*i+{j}]"


if s < 0:
    print(typestr, f"t[{tcnt}];")
print(f"for(i = 0;i < N / {tcnt};++i)" "{")

if s < 0:
    for j in range(tcnt):
        print(f"{getj(j, True)} = {ename} - {getj(j, False)};")


j = 0
jbit = 0
k = 0
kbit = 0
while j < tcnt and k < rcnt:
    barornot = "|" if kbit else " "
    sft = f"<< {kbit - jbit}" if jbit < kbit else f">> {jbit - kbit}"
    print(f"r[{rcnt}*i+{k}] {barornot}= {getj(j, s < 0)} {sft};")
    diff = min(blen - jbit, 8 - kbit)
    jbit += diff
    kbit += diff
    if jbit == blen:
        j += 1
        jbit = 0
    if kbit == 8:
        k += 1
        kbit = 0


print("}")

print("unpack part")

print(f"for(i = 0;i < N / {tcnt};++i)" "{")
j = 0
jbit = 0
k = 0
kbit = 0
while j < tcnt and k < rcnt:
    barornot = "|" if jbit else " "
    sft = f"<< {jbit - kbit}" if kbit < jbit else f">> {kbit - jbit}"
    typeornot = f"({typestr})" if kbit < jbit else ""
    print(f"r->coeffs[{tcnt}*i+{j}] {barornot}= {typeornot}a[{rcnt}*i+{k}] {sft};")
    diff = min(blen - jbit, 8 - kbit)
    jbit += diff
    kbit += diff
    if jbit == blen:
        print(f"r->coeffs[{tcnt}*i+{j}] &= {andstr};")
        print()
        j += 1
        jbit = 0
    if kbit == 8:
        k += 1
        kbit = 0
if s < 0:
    for j in range(tcnt):
        print(f"r->coeffs[{tcnt}*i+{j}] = {ename} - r->coeffs[{tcnt}*i+{j}];")


print("}")
