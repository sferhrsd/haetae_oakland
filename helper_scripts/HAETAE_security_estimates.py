from sub_procedures.MSIS_security import MSIS_summarize_attacks, MSISParameterSet
from sub_procedures.MLWE_security import MLWE_summarize_attacks, MLWEParameterSet
from math import sqrt, log, pi, floor, ceil
from scipy.special import betaincinv, gammaln
from sub_procedures.entropy_coordinate_hyperball import compute_entropy_coordinate, compute_f_table
from sub_procedures.tau_computer import tau_computer
from sub_procedures.gamma_Estimator import gamma_estimate
import sys

# This script is a modification of the one for Dilithium, available at
# https://github.com/pq-crystals/security-estimates
# Please run python3 HAETAE-estimates.py --help for a list of commands.

#Select which estimators should be ran
#The key recovery estimator is much slower than the others and should be toggled off when k,l and eta are unchanged
weak_uf = True
strong_uf = True
key_recovery = True
size = True
fast_LWE = False

if "--no-weak_uf" in sys.argv:
    weak_uf = False
if "--no-strong_uf" in sys.argv:
    strong_uf = False
if "--no-key_recovery" in sys.argv:
    key_recovery = False
if "--no-size" in sys.argv:
    size = False
if "--fast_LWE" in sys.argv or "-f" in sys.argv:
    fast_LWE = True
if len(sys.argv)<=1 or "--help" in sys.argv:
    print("Security Estimator for HAETAE signature.\nOPTIONS:\n--record_xxx: Adds the record parameter sets for xxx security (xxx can be 120, 180 or 260) to the list of parameters.\n--param=\"n=256 q=122889 k=2 l=2 eta=3 security=120 adapt=1\": adds the following parameter set to the list. Multiple parameter sets can be input this way. security must be 120, 180 or 260. If adapt is set to 1, k, l and eta may be increased until target security is reached.\n--no-weak_uf: skips the weak unforgeability hardness computation.\n--no-strong_uf: skips the strong unforgeability hardness computation.\n--no-key_recovery: skips the key recovery hardness computation.\n--no-size: skips the expected verification key and signature sizes computation.\n--fast_LWE or -f: skips the dual attack for LWE parameters. Makes the script a lot faster. Note that this may overestimate the LWE cost by a few bits.")
    exit()

###########################################################################
######################## PRELIMINARIES ####################################
###########################################################################


#This defines, for each security level, the necessary entropy in the hash function.
entropy = { 120 : 192, 180 : 225, 260 : 257 }


class HAETAEParameterSet(object):
    def __init__(self, q, n=256, k=1, l=1, gamma = 0, eta=5, tau=0, security=120, alpha = 512, M = 4, keygen_rate = 50, d = -1):
        """
        Defines a Haetae parameter set
        """
        # If tau = 0, it will be computed using the "entropy" dictionary above.

        #Entropy of hash function
        self.tau = tau
        if tau == 0:
            self.tau = tau_computer(n, entropy[security])
            print("Determining best value for tau: "+str(self.tau)+".\n")
        
        #KeyGen parameters: truncation, rejection rate and gamma bound
        self.d = d
        self.keygen_rate = keygen_rate
        self.gamma = gamma
        if gamma==0:
            print("Gamma not given as input. Computing a rough estimation of it...")
            self.gamma = gamma_estimate(10000, n, k+l, eta, self.tau, [keygen_rate], d, k)[0]
            print("Done! Found gamma="+str(self.gamma))

        #Ring dimension and modulus
        self.n   = n
        self.q   = q

        #Expected number of iterations (must be consistent with N_lower_bound.py)
        M=M/1.001
        self.M = M

        #LWE parameters
        self.l   = l
        self.eta = eta
        self.k = k
        self.beta = self.gamma*sqrt(self.tau)
        self.m = self.n*(self.k+self.l)
        self.gamma1 = self.beta/sqrt((M/2)**(2/self.m)-1)
        self.alpha = alpha
        self.B = self.gamma1+sqrt(self.k*self.n)*(self.alpha/4+1)+sqrt(self.m)/2 

        #Hint truncation parameter
        self.alpha = alpha
        #if d<0:
        #    self.gamma2 = 0
        #else:
        #    #self.gamma2 = tau*2**(d-1)
        #    self.gamma2 = int(sqrt(2*log(2*self.tau/0.01)*self.tau/12)*2**(d-1))+1
        # SIS ell_2 bound for unforgeability, using selftargetMSIS
        self.zeta = self.B
        #if d>=0:
        #    self.zeta += sqrt(self.n*self.k)*(2**(d-1))
        # SIS ell_2 bound for strong unforgeability
        self.zeta_prime = 2*self.B
        if(self.zeta_prime>self.q):
            print("Warning: SIS bound too big for strong unforgeability.")



# Schemes parameters
####################
n= 256

all_params = []

#Parse input parameters if given
for s in sys.argv:
    if s[:len("--param=")] == "--param=":
        args = s[len("--param="):].split(',') #Remove the beginning and final "
        print("Currently parsing",args,"\n")
        values = {"q" : "64513", "n" : "256", "k" : "2", "l" : "4", "security" : "120", "eta" : "1", "M" : "5.5", "keygen_rate" : "25", "d" : "4", "alpha" : "256", "gamma" : "0"}
        for opt in args:
            print(opt)
            values[opt.split("=")[0]] = opt.split("=")[1]
            print(values[opt.split("=")[0]])
            print(opt.split("=")[0])
        print(values["l"])
        all_params += [("HAETAE "+str(s), HAETAEParameterSet(int(values["q"]), int(values["n"]), k=int(values["k"]), l=int(values["l"]), eta=int(values["eta"]), security=int(values["security"]), M=float(values["M"]), keygen_rate = int(values["keygen_rate"]), d = int(values["d"]), alpha=int(values["alpha"]), gamma = float(values["gamma"])))]


######################
# Hardcoded parameters
######################

if "--15_bits" in sys.argv:
    all_params += [("HAETAE 15 bits", HAETAEParameterSet(32257, n, k=2, l=4, eta=1, security = 120, M = 6, keygen_rate = 10, d = -1, alpha = 128, gamma = 42.35))]
if "--record_120" in sys.argv:
#Current Implem
    all_params += [("HAETAE Medium", HAETAEParameterSet(64513, n, k=2, l=4, eta=1, security = 120, M = 6, keygen_rate = 10, d = 1, alpha = 512, gamma = 46.59))]

# Params with various truncation removed to decrease M
#    all_params += [("HAETAE Medium", HAETAEParameterSet(64513, n, k=2, l=4, eta=1, security = 120, M = 4.8, keygen_rate = 15, d = 1, alpha = 256, gamma = 48))]
#    all_params += [("HAETAE Medium", HAETAEParameterSet(64513, n, k=2, l=4, eta=1, security = 120, M = 5, keygen_rate = 25, d = 0, alpha = 512, gamma = 43.5))]
#    all_params += [("HAETAE Medium", HAETAEParameterSet(64513, n, k=2, l=4, eta=1, security = 120, M = 4, keygen_rate = 25, d = 0, alpha = 256, gamma = 43.5))]

#    all_params += [("HAETAE Medium", HAETAEParameterSet(64513, n, k=3, l=4, eta=1, security = 120, M = 4, keygen_rate = 25, d = 2, alpha = 1024, gamma = 70.58))]

if "--record_180" in sys.argv:
    all_params += [("HAETAE Recommended", HAETAEParameterSet(64513, n, k=3, l=6, eta=1, security = 180, M = 5, keygen_rate = 25, d = 1, alpha = 512, gamma = 56))]
#    all_params += [("HAETAE Recommended", HAETAEParameterSet(64513, n, k=4, l=5, eta=2, security = 180, M = 5, keygen_rate = 10, d = 1, alpha = 512, gamma = 90.81))]
#Record for 260609
#    all_params += [("HAETAE Recommended", HAETAEParameterSet(260609, n, k=3, l=6, eta=1, security = 180, M = 4.5, keygen_rate = 25, d = 2, alpha = 2048, gamma = 72.43))]
if "--record_260" in sys.argv:
    all_params += [("HAETAE Very High", HAETAEParameterSet(64513, n, k=4, l=7, eta=1, tau=128, security=260, M = 6, keygen_rate = 10, alpha = 256, d = -1, gamma = 55.13))]
#Record for 260609
#    all_params += [("HAETAE Very High", HAETAEParameterSet(260609, n, k=4, l=7, eta=1, tau=128, security=260, M = 5, keygen_rate = 10, alpha = 2048, d = 1, gamma = 60.75))]


#########################
# Conversion to MSIS/MLWE
#########################

def HAETAE_to_MSIS(dps, strong_uf = False):
    if strong_uf:
        return MSISParameterSet(dps.n, dps.k + dps.l, dps.k, dps.zeta_prime, dps.q, norm="l2")
    return MSISParameterSet(dps.n, dps.k + dps.l, dps.k, dps.zeta, dps.q, norm="l2")


def HAETAE_to_MLWE(dps):
    #Also valid for Ring case as solving the NTRU instance is solving hf+g=0, an LWE instance
    return MLWEParameterSet(dps.n, max(dps.l-1,1), dps.k, dps.eta, dps.q, distr="uniform")

text_SIS = ["BKZ block-size $b$ to break SIS","Best Known Classical bit-cost","Best Known Quantum bit-cost","Best Plausible bit-cost"]
text_LWE = ["BKZ block-size $b$ to break LWE","Best Known Classical bit-cost","Best Known Quantum bit-cost","Best Plausible bit-cost"]


##################
# Size Computation
##################

def HAETAE_Entropy(dps):
    """
    Computes the optimal expected size of a signature by replacing the encoding of z by its entropy
    Underestimate by a few bytes the real value, even more when alpha grows bigger.
    """
    #Size of tilde{c} is n bits
    size_c = dps.n 
    m = dps.n*(dps.l+dps.k)
    (p, _) = compute_f_table(m, floor(dps.gamma1), int(log(dps.alpha,2)))
    hint_entropy = sum([-(l+1e-16)*log(l+1e-16) for l in p])
    #print(hint_entropy*dps.n*dps.k/8)
    size_c_h = size_c + hint_entropy*dps.n*dps.k #Number of bits necessary to represent the hint, i.e. ~High(z_2)
    return size_c_h +(compute_entropy_coordinate(m,floor(dps.gamma1)))*(dps.n*dps.l) #Number of bits necessary to represent z_1, assuming entropic coding

def HAETAE_PK_Size(dps):
    """
    Computes the expected size of a verification key. This does not depend on the distribution used.
    """
    if dps.d <0:
        return (256 + dps.k*dps.n*(int(log(dps.q,2)+1))) #The public vector is computed mod q
    return (256 + dps.k*dps.n*(int(log(dps.q,2)+1)-dps.d)) #The public vector is computed mod q

# The remainder of this script is formatting
#############################################
######################### ANALYSIS AND REPORT 
#############################################


table_weak_SIS   = [len(all_params)*[0] for i in range(4)]
table_strong_SIS = [len(all_params)*[0] for i in range(4)]
table_LWE        = [len(all_params)*[0] for i in range(4)]
table_size       = [0 for i in range(len(all_params))]
table_entropy    = [0 for i in range(len(all_params))]
table_pk         = [0 for i in range(len(all_params))]


#For each selected scheme, build the estimate cost of selected attacks
j = 0
for (scheme, param) in all_params:
    print("\n"+scheme)
    print(param.__dict__)
    print("")
    if weak_uf:
        print("=== WEAK UF")
        v = MSIS_summarize_attacks(HAETAE_to_MSIS(param))
        for i in range(4):
            table_weak_SIS[i][j] = v[i]
    if strong_uf:
        print("=== STRONG UF")
        v = MSIS_summarize_attacks(HAETAE_to_MSIS(param, strong_uf=True))
        for i in range(4):
            table_strong_SIS[i][j] = v[i]
    if key_recovery:
        print("=== SECRET KEY RECOVERY")
        v = MLWE_summarize_attacks(HAETAE_to_MLWE(param),fast_LWE)
        for i in range(4):
            table_LWE[i][j] = v[i]
    if size:
        print("=== SIGNATURE SIZE")
        #table_size[distribution][j] = HAETAE_Signature_Size(param)
        #print(table_size[distribution][j])
        table_entropy[j] = HAETAE_Entropy(param)
        print(table_entropy[j]/8)
        print("=== PK SIZE")
        table_pk[j] = HAETAE_PK_Size(param)
        print(table_pk[j]/8)
    j+=1

#Print the generated LaTeX table
print("HAETAE TABLE")
print("========================")
print("\\hline")
print("$q$"+"".join([" & "+str(dps[1].q) for dps in all_params]))
print("\\\\")
print("$M$"+"".join([" & "+str(dps[1].M*1.001) for dps in all_params]))
print("\\\\")
print("$\\gamma$"+"".join([" & "+str(dps[1].gamma) for dps in all_params]))
print("\\\\")
print("KeyGen Acceptance Rate"+"".join([" & "+str(dps[1].keygen_rate/100) for dps in all_params]))
print("\\\\")
print("$\\beta=S\\sqrt{\\tau}$"+"".join([" & "+"{:5.2f}".format(dps[1].beta) for dps in all_params]))
print("\\\\")
print("$B$"+"".join([" & "+"{:6.2f}".format(sqrt(dps[1].gamma1**2+dps[1].beta**2)) for dps in all_params]))
print("\\\\")
print("$B'$"+"".join([" & "+"{:6.2f}".format(dps[1].gamma1) for dps in all_params]))
print("\\\\")
print("$B''$"+"".join([" & "+"{:6.2f}".format(dps[1].zeta) for dps in all_params]))
print("\\\\")
print("$(k,\\ell)$"+"".join([" & ("+str(dps[1].k)+","+str(dps[1].l)+")" for dps in all_params]))
print("\\\\")
print("$\\eta$"+"".join([" & "+str(dps[1].eta) for dps in all_params]))
print("\\\\")
print("$\\tau$"+"".join([" & "+str(dps[1].tau) for dps in all_params]))
print("\\\\")
print("$\\alpha$"+"".join([" & "+str(dps[1].alpha) for dps in all_params]))
print("\\\\")
print("$d$"+"".join([" & "+str(dps[1].d) for dps in all_params]))
print("\\\\")
print("\\hline")
for j in range(4):
    print(text_SIS[j]+"".join([" & "+str(table_weak_SIS[j][i])+" ("+str(table_strong_SIS[j][i])+")" for i in range(len(all_params))]))
    print("\\\\")
print("\\hline")
for j in range(4):
    print(text_LWE[j]+"".join([" & "+str(table_LWE[j][i]) for i in range(len(all_params))]))
    print("\\\\")
print("\\hline")
print("Signature size with rANS"+"".join([" & "+str(int(table_entropy[i]/8)) for i in range(len(all_params))]))
print("\\\\")
print("Public key size"+"".join([" & "+str(int(table_pk[i]/8)) for i in range(len(all_params))]))
print("\\\\")
print("Sum"+"".join([" & "+str(int((table_pk[i]+table_entropy[i])/8)) for i in range(len(all_params))]))
print("\\\\")
print("\\hline")
print("========================")
