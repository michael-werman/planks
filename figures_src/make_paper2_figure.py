# Illustration for Paper II (floating point; not used in any proof).
# Left: a convex body (ellipse), n random points, the spanned lines, and the point of the body
#       farthest from all of them (found on a grid), with the empty disk of radius rho.
# Right: a point z of the body and the chords through it; h_Omega(z) is proportional to the
#        integral of the cubed chord lengths (Corollary 2.2).
import numpy as np, itertools, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
rng=np.random.default_rng(5)
A,B=1.4,1.0                      # ellipse semi-axes
def inside(P): return (P[...,0]/A)**2+(P[...,1]/B)**2<=1
def sample(n):
    out=[]
    while len(out)<n:
        p=rng.uniform([-A,-B],[A,B])
        if inside(p): out.append(p)
    return np.array(out)
def chord(p,u):                  # intersection of line p+tu with the ellipse
    a=(u[0]/A)**2+(u[1]/B)**2; b=2*(p[0]*u[0]/A**2+p[1]*u[1]/B**2); c=(p[0]/A)**2+(p[1]/B)**2-1
    disc=b*b-4*a*c
    if disc<=0: return None
    r=np.sqrt(disc); t1,t2=(-b-r)/(2*a),(-b+r)/(2*a)
    return p+t1*u,p+t2*u
X=sample(7)
g=np.linspace(-A,A,701); G=np.array(np.meshgrid(g,np.linspace(-B,B,501))).reshape(2,-1).T; G=G[inside(G)]
best=np.zeros(len(G))+np.inf
segs=[]
for i,j in itertools.combinations(range(len(X)),2):
    u=X[j]-X[i]; u/=np.linalg.norm(u); nrm=np.array([-u[1],u[0]])
    best=np.minimum(best,np.abs((G-X[i])@nrm)); segs.append(chord(X[i],u))
k=np.argmax(best); z,rho=G[k],best[k]
th=np.linspace(0,2*np.pi,400)
fig,ax=plt.subplots(1,2,figsize=(9.5,3.9))
for a in ax:
    a.plot(A*np.cos(th),B*np.sin(th),color="black",lw=1.2); a.set_aspect("equal"); a.axis("off")
for s in segs: ax[0].plot([s[0][0],s[1][0]],[s[0][1],s[1][1]],color="#4a6fa5",lw=0.8)
ax[0].plot(X[:,0],X[:,1],"o",color="black",ms=3.5)
ax[0].add_patch(plt.Circle(z,rho,color="#e07b39",alpha=0.35)); ax[0].plot(*z,"o",color="#e07b39",ms=3.5)
ax[0].set_title(r"spanned lines and the largest empty disk ($\rho$)",fontsize=9)
z2=np.array([0.55,-0.35])
for t in np.linspace(0,np.pi,12,endpoint=False):
    s=chord(z2,np.array([np.cos(t),np.sin(t)]))
    ax[1].plot([s[0][0],s[1][0]],[s[0][1],s[1][1]],color="#4a6fa5",lw=0.8)
ax[1].plot(*z2,"o",color="#e07b39",ms=4); ax[1].annotate(r"$z$",z2+np.array([0.05,0.06]),fontsize=10)
ax[1].set_title(r"chords through $z$: $h_\Omega(z)\propto\int_0^\pi L_z(\theta)^3\,d\theta$",fontsize=9)
plt.tight_layout(rect=(0,0,1,0.93)); plt.savefig("../paper2/figures/fig_lines.png",dpi=200)
print("rho (grid, floating point) =",rho)
