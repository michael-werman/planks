# EXACT (sympy): limit constant for j-flats spanned by j+1 points with density f_beta ~ (1-|x|^2)^beta in the unit
# ball of R^d (beta = 0: uniform), from c = (j+1)! d (j+1)/2 * B(a/2,k+1) B(a/2,j/2) / B(a/2,k+j/2),
# a = d-j, k = j(d+1)/2 + beta(j+1); plus a floating-point Monte Carlo check of h at |z| = 0 and 1.
import sympy as sp
def c(d,j,beta=0):
    a=sp.Integer(d-j); k=sp.Rational(j*(d+1),2)+beta*(j+1)
    B=lambda x,y: sp.gamma(x)*sp.gamma(y)/sp.gamma(x+y)
    return sp.nsimplify(sp.simplify(sp.factorial(j+1)*d*(j+1)/2*B(a/2,k+1)*B(a/2,sp.Rational(j,2))/B(a/2,k+sp.Rational(j,2))))
def h(d,j,r,beta=0):  # h at |z| = r (floating point quadrature over the Beta law)
    import mpmath as mp
    a=d-j; k=mp.mpf(j*(d+1))/2+beta*(j+1)
    Bf=lambda x,y: mp.beta(x,y)
    dens=lambda x: x**(a/2-1)*(1-x)**(j/2-1)/Bf(a/2,j/2)
    return 2/(a*Bf(a/2,k+1))*mp.quad(lambda x:(1-r*r*x)**k*dens(x),[0,1])
if __name__=="__main__":
    for d in range(2,6):
        print(d,[ (j, c(d,j), float(c(d,j))) for j in range(1,d)])
    import numpy as np
    rng=np.random.default_rng(3)
    def ball(n,d):
        g=rng.standard_normal((n,d)); g/=np.linalg.norm(g,axis=1)[:,None]; return g*rng.uniform(0,1,(n,1))**(1/d)
    for d,j in [(3,1),(4,1),(4,2)]:
        a=d-j; N=2_000_000; t=0.02
        P=[ball(N,d) for _ in range(j+1)]
        Q=np.stack([p-P[0] for p in P[1:]],2)            # N x d x j
        Qo,_=np.linalg.qr(Q)                               # orthonormal basis of direction
        for r in (0.0,1.0):
            z=np.zeros(d); z[0]=r; w=z-P[0]
            perp=w-np.einsum('nij,nj->ni',Qo,np.einsum('nij,ni->nj',Qo,w))
            est=(np.linalg.norm(perp,axis=1)<t).mean()/t**a
            print(f"d={d} j={j} |z|={r}: MC {est:.3f}, formula {float(h(d,j,r)):.3f}")
