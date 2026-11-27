# FLOATING POINT exploration (not a proof): limit constant c_beta for the beta density
# f = const*(1-|x|^2)^beta on the unit ball of R^d (beta >= 0; beta = 0 is uniform), from the formula
# c = (d^2 d!/2) J((d-3)/2) J(k)/J(k+(d-3)/2), k = (d^2-1)/2 + d*beta, J(a) = sqrt(pi) Gamma(a+1)/Gamma(a+3/2),
# plus a Monte Carlo check of h_f at |z| = 1 for d = 2.
import numpy as np
from math import gamma, pi, sqrt, factorial
J=lambda a: sqrt(pi)*gamma(a+1)/gamma(a+1.5)
def c(d,b):
    k=(d*d-1)/2+d*b
    return d*d*factorial(d)/2*J((d-3)/2)*J(k)/J(k+(d-3)/2)
def m2(b):  # min h_f for d=2
    k=1.5+2*b; return 2*J(k-0.5)/(J(-0.5)*J(k))
if __name__=="__main__":
    for d in (2,3,4):
        print(d,[round(c(d,b),4) for b in (0,0.25,0.5,1,2,5,20)])
    rng=np.random.default_rng(1)
    def sample(n,b):  # radial: |x|^2 ~ Beta(1, b+1) in d=2
        r=np.sqrt(rng.beta(1,b+1,n)); t=rng.uniform(0,2*np.pi,n)
        return np.stack([r*np.cos(t),r*np.sin(t)],1)
    for b in (0,1,2):
        N=4_000_000; t=0.004; P=sample(N,b); Q=sample(N,b); z=np.array([1.0,0])
        u=Q-P; nrm=np.stack([-u[:,1],u[:,0]],1)/np.linalg.norm(u,axis=1)[:,None]
        dist=np.abs(((z-P)*nrm).sum(1)); est=(dist<t).mean()/t
        print(f"beta={b}: MC h(1,0)={est:.4f} +- {np.sqrt((dist<t).mean()/N)/t:.4f}, formula {m2(b):.4f}")
