# FLOATING POINT exploration (not a proof): h(z) = (2/(3|K|^2)) * int_0^pi L_z(theta)^3 dtheta for convex polygons,
# minimized over a grid; reports the minimizer and c = 8/min h.
import numpy as np, sys
def chord(P,z,th):
    d=np.array([np.cos(th),np.sin(th)]); lo,hi=-np.inf,np.inf
    m=len(P)
    for i in range(m):
        a,b=P[i],P[(i+1)%m]; nrm=np.array([b[1]-a[1],a[0]-b[0]])  # outward for CCW
        den=nrm@d; num=nrm@(a-z)
        if abs(den)<1e-15:
            if num<-1e-15: return 0.0
            continue
        t=num/den
        if den>0: hi=min(hi,t)
        else: lo=max(lo,t)
    return max(hi-lo,0.0)
def h(P,z,K=4000):
    th=(np.arange(K)+0.5)*np.pi/K
    area=0.5*abs(sum(P[i][0]*P[(i+1)%len(P)][1]-P[(i+1)%len(P)][0]*P[i][1] for i in range(len(P))))
    return 2/(3*area**2)*np.mean([chord(P,z,t)**3 for t in th])*np.pi
shapes={'square':np.array([[0,0],[1,0],[1,1],[0,1.]]),
        'triangle':np.array([[0,0],[1,0],[0.5,np.sqrt(3)/2]])}
if __name__=="__main__":
  for name,P in shapes.items():
      best=(9,None); xs=np.linspace(0,1,41)
      for x in xs:
          for y in np.linspace(0,1,41):
              z=np.array([x,y])
              if all((np.array([P[(i+1)%len(P)][1]-P[i][1],P[i][0]-P[(i+1)%len(P)][0]])@(z-P[i]))<=1e-12 for i in range(len(P))):
                  v=h(P,z,800); best=min(best,(v,tuple(z)))
      v=h(P,np.array(best[1]),20000)
      print(name,'min h ~',round(v,5),'at',best[1],' c=8/min h ~',round(8/v,4),' centroid h ~',round(h(P,P.mean(0),4000),4))
