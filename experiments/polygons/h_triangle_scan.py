# FLOATING POINT exploration (not a proof): minimize h(z) = (2/(3|K|^2)) int_0^pi L_z^3 over a triangle,
# scanning the boundary finely and the interior on a grid. Uses chord/h from h_polygon.py.
import numpy as np
from h_polygon import h
def scan(P,nb=200,ng=40,K=600):
    best=(9,None,'')
    for i in range(3):
        a,b=P[i],P[(i+1)%3]
        for s in np.linspace(0,1,nb+1):
            z=a+s*(b-a); v=h(P,z,K); best=min(best,(v,tuple(np.round(z,4)),'edge %d s=%.3f'%(i,s)))
    bi=9
    for u in np.linspace(0,1,ng+1)[1:-1]:
        for w in np.linspace(0,1,ng+1)[1:-1]:
            if u+w<1:
                z=P[0]+u*(P[1]-P[0])+w*(P[2]-P[0]); bi=min(bi,h(P,z,K))
    v=h(P,np.array(best[1]),20000)
    print('min on boundary h~%.5f at %s (%s); c=8/h~%.4f; min over interior grid ~%.4f'%(v,best[1],best[2],8/v,bi))
scan(np.array([[0,0],[1,0],[0.5,np.sqrt(3)/2]]))
scan(np.array([[0,0],[1,0],[0.05,0.3]]))
