import numpy as np
import scipy.fft as dft
from math import sqrt

def gamma_estimate(N, n, dim, eta, tau, rate, d, k):
    """
    Samples N secret keys for random uniform matrices.
    Compute the gamma value, as given in the documentation.
    returns all the rate/100-th quantile of these values.
    """
    res = []
    i_max = n//tau
    leftover = n % tau

    identity = np.array([1]+[0 for i in range(2*n-1)])
    for loop in range(N):
        #Generate a secret
        s1 = [np.concatenate((np.random.randint(-eta,eta+1,size=n), np.array([0 for i in range(n)]))) for i in range(dim)]
        s1[0] = identity
        #Generate a public key (assumed uniform here)
        if d>1:
            #Do the High/Low bits truncation
            b0 = [np.concatenate(([np.random.randint(-2**(d-1)+1,2**(d-1)+1) for j in range(n)], np.array([0 for j in range(n)]))) for i in range(k)]
            #If we have the extremal value, it has probability one half of being flipped
            b0 = [np.array([val if val!= 2**(d-1) else (-1)**np.random.randint(0,2)*val for val in b0[i]]) for i in range(k)]
            s1 = [s1[i] if i<dim-k else s1[i]-b0[i-dim+k] for i in range(dim)]
        elif d==1:
            b0 = [np.concatenate(([2*np.random.randint(0,2)-1 if np.random.randint(0,2)==0 else np.random.randint(-2**(d-1)+1,2**(d-1)) for i in range(n)], np.array([0 for i in range(n)]))) for i in range(k)]
            s1 = [s1[i] if i<dim-k else s1[i]-b0[i-dim+k] for i in range(dim)]
        #Compute the canonical embeddings and their norm
        y = [dft.fft(s1[i])[1::2] for i in range(len(s1))]
        norm_y = [np.linalg.norm([abs(y[i][j]) for i in range(len(y))]) for j in range(len(y[0]))]
        largest = max(norm_y)
        #We now compute the bound
        sorted_y = sorted(norm_y, reverse = True)
        res.append(sqrt(tau**2*sum([x**2 for x in sorted_y[:i_max]])+(leftover*tau)*sorted_y[i_max]**2)/sqrt(n*tau))
    return([np.nanquantile(res, r/100) for r in rate])

#print("Best 25%, 50% and max values for the two bounds: "+str(S_estimate(10000, 256, 6, 1, 58, [10,25,50], 1, 2)))
