#Input the ring dimension and target entropy.
#Output the smallest tau such that c has sufficient min entropy.

import math

def tau_computer(n,target_entropy):
    tau = 0
    while(math.log(math.comb(n,tau),2))<target_entropy and tau<n:
        tau +=1
    if tau==n:
        print("Warning, overflow.")
        #tau = 0
        #val = math.comb(n,0)
        #while val<2**(256-target_entropy):
        #    tau += 1
        #    val += math.comb(n,tau)
    return tau
