# FLOATING POINT check of h(z) = lim P(dist(z, line U1U2) < t)/t in the unit disc
# against the formula h(z) = (2/(3|Omega|^2)) * int_0^pi L_z(theta)^3 dtheta.
import numpy as np
rng=np.random.default_rng(3)
def disc(n):
    r=np.sqrt(rng.random(n)); a=2*np.pi*rng.random(n); return np.c_[r*np.cos(a),r*np.sin(a)]
def formula(rho):
    th=np.linspace(0,np.pi,200001); L=2*np.sqrt(np.clip(1-rho**2*np.sin(th)**2,0,None))
    return 2/(3*np.pi**2)*np.trapezoid(L**3,th)
N=4*10**7; t=0.002
for rho in [0.0,0.5,0.9,1.0]:
    z=np.array([rho,0.0]); U=disc(N); V=disc(N); d=V-U
    dist=np.abs(d[:,0]*(z[1]-U[:,1])-d[:,1]*(z[0]-U[:,0]))/np.linalg.norm(d,axis=1)
    print(rho, (dist<t).mean()/t, formula(rho))
print("limit constant 8/h(boundary) =", 8/formula(1.0), "9pi^2/8 =", 9*np.pi**2/8)
