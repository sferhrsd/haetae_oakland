from gamma_Estimator import gamma_estimate

N   = 10000
n   = 256
dim = 7
eta = 1
tau = 58
#d   = 3
k   = 3

print(gamma_estimate(N, n, dim, eta, tau, [10 ,25,50], 0, k))
print(gamma_estimate(N, n, dim, eta, tau, [1, 5, 10, 25, 50], 1, k))
print(gamma_estimate(N, n, dim, eta, tau, [10 ,25,50], 2, k))
print(gamma_estimate(N, n, dim, eta, tau, [10 ,25,50], 3, k))
print(gamma_estimate(N, n, dim, eta, tau, [10 ,25,50], 4, k))
print(gamma_estimate(N, n, dim, eta, tau, [10 ,25,50], 5, k))
