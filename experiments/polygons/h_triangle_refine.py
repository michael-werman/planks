# FLOATING POINT exploration (not a proof): refines the boundary minimum of h for the lopsided triangle
# (0,0),(1,0),(0.05,0.3), comparing the bottom edge with the edge from (1,0) to (0.05,0.3).
# h(z) = (2/(3|K|^2)) int_0^pi L_z(theta)^3 dtheta, midpoint rule with K theta nodes (vectorized), golden-section search.
import numpy as np
P=np.array([[0,0],[1,0],[0.05,0.3]])
u,v=P[1]-P[0],P[2]-P[0]; area=0.5*abs(u[0]*v[1]-u[1]*v[0])
def h(z,K):
    th=(np.arange(K)+0.5)*np.pi/K; d=np.stack([np.cos(th),np.sin(th)],1)
    lo=np.full(K,-np.inf); hi=np.full(K,np.inf)
    for i in range(3):
        a,b=P[i],P[(i+1)%3]; nrm=np.array([b[1]-a[1],a[0]-b[0]])
        den=d@nrm; num=nrm@(a-z)
        with np.errstate(divide='ignore',invalid='ignore'): t=num/den
        hi=np.where(den>1e-15,np.minimum(hi,t),hi); lo=np.where(den<-1e-15,np.maximum(lo,t),lo)
    L=np.clip(hi-lo,0,None)
    return 2/(3*area**2)*np.mean(L**3)*np.pi
def gs(f,a,b,it=60):
    g=(np.sqrt(5)-1)/2; c=b-g*(b-a); e=a+g*(b-a); fc=f(c); fe=f(e)
    for _ in range(it):
        if fc<fe: b,e,fe=e,c,fc; c=b-g*(b-a); fc=f(c)
        else: a,c,fc=c,e,fe; e=a+g*(b-a); fe=f(e)
    return (a+b)/2
if __name__=="__main__":
    for K in (2000,200000,1000000):
        x=gs(lambda x:h(np.array([x,0.]),K),0.2,0.5); hb=h(np.array([x,0.]),K)
        e=lambda s:P[1]+s*(P[2]-P[1])
        s=gs(lambda s:h(e(s),K),0.5,0.8); he=h(e(s),K)
        print(f"K={K}: bottom edge min h={hb:.6f} at ({x:.4f},0), c={8/hb:.5f};  slanted edge min h={he:.6f} at {np.round(e(s),4)}, rel gap={(he-hb)/hb:.2e}")
