q = 12289 # prime for ring
N = 256 # polynomial length
R = 2 ** 32 # montgomery factor
ROOT_OF_UNITY = 3

def modpm(num): # rounding to (-q / 2, q / 2]
    rd = num % q
    if rd > q / 2:
        return -(q-rd)
    else:
        return rd





l = N // 2 #start length(N / 2)
lis = [0]
nowarr = [128]
nxtarr = []
while l > 0:
    for i in nowarr:
        lis.append(modpm((ROOT_OF_UNITY ** i) * R))
        nxtarr.append(i // 2)
        nxtarr.append(i // 2 + 128)
    nowarr = nxtarr
    nxtarr = []
    l >>= 1
print("{", end = "")
print(*lis, sep = ', ', end = "")
print("};")
