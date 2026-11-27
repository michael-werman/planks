# FLOATING POINT. rho for n uniform points in unit disc, lines through pairs; grid lower bound + Lipschitz upper bound.
import numpy as np, sys, itertools
rng=np.random.default_rng(11)
def disc(m):
    a=rng.uniform(0,2*np.pi,m); r=np.sqrt(rng.random(m)); return np.c_[r*np.cos(a),r*np.sin(a)]
def rho(X,h):
    g=np.arange(-1,1+h,h); G=np.array(np.meshgrid(g,g)).reshape(2,-1).T; G=G[(G**2).sum(1)<=1]
    th=np.arange(0,2*np.pi,h/2); G=np.r_[G,np.c_[np.cos(th),np.sin(th)]]
    i,j=np.array(list(itertools.combinations(range(len(X)),2))).T
    e=X[j]-X[i]; nn=np.c_[-e[:,1],e[:,0]]; nn/=np.linalg.norm(nn,axis=1)[:,None]; c=(nn*X[i]).sum(1)
    best=0
    for k in range(0,len(G),20000):
        d=np.abs(G[k:k+20000]@nn.T-c).min(1); best=max(best,d.max())
    return best
for n,h,T in [(8,0.002,40),(12,0.002,40),(20,0.002,30),(30,0.0015,20),(40,0.0015,20)]:
    lo=[];
    for _ in range(T):
        lo.append(rho(disc(n),h))
    lo=np.array(lo); s=n*n/np.log(n)
    print('n=%d median n^2 rho/log n in [%.2f, %.2f]  (quartiles %.2f %.2f)'%(n,np.median(lo)*s,(np.median(lo)+h)*s,np.quantile(lo,.25)*s,np.quantile(lo,.75)*s),flush=True)
