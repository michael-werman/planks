# EXACT (sympy) and floating-point exploration: limit constant for the radial densities f_a = (1+a|x|^2)/(pi(1+a/2))
# on the unit disk (d=2). g(p) = chord weight int int f f |s-t| over the chord at offset p; b_2 = 1/2;
# h(z) = int_0^{2pi} g(|z| cos th) d th (Lemma 3.3(c)); c = 8/min h.
import sympy as sp
a,p,s,t,th,r=sp.symbols('a p s t theta r',real=True)
w=sp.sqrt(1-p**2)
Z=sp.pi*(1+a/2)
F=lambda u:(1+a*(p**2+u**2))/Z
# |s-t| integral over the square: 2 * integral over t<s
g=2*sp.integrate(sp.integrate(F(s)*F(t)*(s-t),(t,-w,s)),(s,-w,w))
g=sp.simplify(g)
if __name__=="__main__":
    print('g(p)=',sp.factor(g))
    hb=sp.simplify(2*sp.integrate(g.subs(p,sp.cos(th)).subs(sp.sqrt(1-sp.cos(th)**2),sp.sin(th)),(th,0,sp.pi)))
    hb=sp.simplify(sp.integrate(sp.expand(sp.simplify(g.subs(p,sp.cos(th)))),(th,0,2*sp.pi)))
    print('h on the circle =',sp.simplify(hb))
    c=sp.simplify(8/hb); print('c(a)=',c, ' c(0)=',sp.simplify(c.subs(a,0)))
    print("dc/da at 0 =",sp.simplify(sp.diff(c,a).subs(a,0)), float(sp.diff(c,a).subs(a,0)))
    for av in [0,0.25,0.5,1,2,5,20]:
        print(av, float(c.subs(a,av)))

def mc_check(N=4_000_000,t=0.004,seed=2):
    # FLOATING POINT Monte Carlo of t^{-1} Pr(dist(z,H)<t) at |z|=1 for f = (1+2|x|^2)/(2 pi); formula 27312/(3675 pi^2).
    import numpy as np
    rng=np.random.default_rng(seed)
    def sample(n):
        u=rng.uniform(0,1,n); R=np.sqrt((-1+np.sqrt(1+8*u))/2); T=rng.uniform(0,2*np.pi,n)
        return np.stack([R*np.cos(T),R*np.sin(T)],1)
    P,Q=sample(N),sample(N); v=Q-P; nrm=np.stack([-v[:,1],v[:,0]],1)/np.linalg.norm(v,axis=1)[:,None]
    hit=(np.abs(((np.array([1.0,0])-P)*nrm).sum(1))<t).mean()
    return hit/t, np.sqrt(hit/N)/t, 27312/(3675*np.pi**2)
