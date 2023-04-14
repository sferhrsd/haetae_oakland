###################################
# This supporting script for HAETAE
# aims at computing a lower bound
# on the discretization parameter N
###################################

from math import sqrt,ceil

# Which values for M0 and c?
M0 = 1.001
c = 1.001

#Our target number of iterations is 4.
M=4/c
print("Value for M",M,"\n")

for (m,t,sec) in [(256*6,293.51,120),(256*9,385,180),(256*11,457.01,260)]:
    r = t/sqrt((M/2)**(2/m)-1)
    rp = sqrt(r**2+t**2)
    
    #See Lemma 2: lower bound for N due to the sampling of y
    N0 = sqrt(m)*(M0**(1/m)+1)/(2*rp*(M0**(1/m)-1))
    #See Lemma 1: lower bound for N due to the constraint on the expected number of iterations of Sign
    N1 = sqrt(m)*(c**(1/m)/r+1/rp)/(2*(c**(1/m)-1))
    
    print("Lower bounds for HAETAE",sec,"\n---------------------------")
    print("Lower bounds:",N0,N1)
    print("Max of the two:",ceil(max(N0,N1)),"\n")
    print("New values for r and r':",r,rp)
