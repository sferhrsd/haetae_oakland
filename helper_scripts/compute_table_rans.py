from sub_procedures.entropy_coordinate_hyperball import compute_f_table
from math import floor, ceil, log

param = [#deg,k, l,     B, a,  n, sec
        (256, 2, 4,  9389, 9, 16, 120), # 120 params
        (256, 3, 6, 17774, 9, 16, 180), # 180 params
        (256, 4, 7, 22344, 8, 16, 260)  # 260 params
        ]

for (deg, k, l, B, a, n, sec) in param:
    print("Computing rANS table for HAETAE",sec)
    m = deg*(k+l)
    (prob, low_prob) = compute_f_table(m,B,a)
    
    #Treatment of low bits: compute their statistical distance with uniform distribution
    dist = sum([abs(i-1/len(low_prob)) for i in low_prob])/2
    print("Statistical distance of low bits with uniform distribution", dist) 
    
    #Treatment of high bits
    prob = [p*2**n for p in prob]
    #Given the probabilities, we multiply them by 2**n.
    #We first check that n is big enough to contain all possible values of the high bits. 
    if(2**n< len(prob)):
        raise Exception("n is too small.")
    #We then round the probabilities to get f. Our heuristic is to set 1 if p<=1 (as we cannot erase values, even if they have extremely low probabilities), and floor everything else.
    f = [1 if p<=1 else floor(p) for p in prob]
    s = sum(f)
    #If s is larger than 2**n, we renormalize and apply the same rounding strategy until we get an admissible table.
    while s>2**n:
        f = [1 if p*2**n/s <=1 else floor(p*2**n/s) for p in f]
        s = sum(f)
    #Now s<=2**n, but we require s=2**n. Add the difference to the largest value
    diff = 2**n - s
    f[f.index(max(f))] += diff
    print(f"Alphabet size = {len(f)}")
    print(f"f = {f}\ncdf = {[sum(f[:i]) for i in range(len(f))]}")
    
    #Finally, we measure the discrepancy between the rounded distribution and the one we computed in order to estimate the encoding loss per coordinate.
    kl = [0 if prob[i]==0 else prob[i]*log(prob[i]/f[i],2)/2**n for i in range(len(prob))]
    print("Average over-cost per coordinate",sum(kl))
    print("")
