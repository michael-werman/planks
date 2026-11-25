# FLOATING POINT exploration for d = 3: per-triple inequality
#   G - (2B/3) e2(x) - (g m0/3) sum(x) - A kappa <= C0 + (1 - b0) m0^2,  x = heights - m0 in [h1-m0, h2-m0]^3
# vertices: j points at h1, 3-j at h2 (j = 0..3). Complete the square in m0. (D, W) = (c cos phi, s sin phi).
import numpy as np, sympy as sp
from scipy.optimize import minimize
h1, h2, m0, B, g = sp.symbols('h1 h2 m0 B g')
G = ((h1+h2)/2 - m0)**2 + ((h2-h1)/2)**2/2
forms = []
for j in range(4):
    x1, x2 = h1 - m0, h2 - m0
    e2 = sp.binomial(j, 2)*x1**2 + j*(3-j)*x1*x2 + sp.binomial(3-j, 2)*x2**2
    sx = j*x1 + (3-j)*x2
    f = sp.expand(G - sp.Rational(2, 3)*B*e2 - g*m0*sx/3)
    f0 = f.subs(m0, 0); L = sp.expand(f.coeff(m0, 1)); k2 = sp.expand(f.coeff(m0, 2))
    forms.append((sp.lambdify((h1, h2, B, g), f0), sp.lambdify((h1, h2, B, g), L), sp.lambdify((B, g), k2)))
print('m0^2 coefficients:', [sp.simplify(sp.expand(G - 0).coeff(m0, 2))])
phi = np.linspace(0, np.pi, 241)[:, None]; psi = np.linspace(0, np.pi/2, 241)[None, :]
c, s = np.cos(psi), np.sin(psi); D = c*np.cos(phi); W = s*np.sin(phi); kap = 1 - c
H1, H2 = D - W, D + W
def worst(A, Bv, gv, b0):
    M = -np.inf
    for f0, L, k2 in forms:
        mu = (1 - b0) - k2(Bv, gv)
        if mu <= 0: return np.inf
        Phi = f0(H1, H2, Bv, gv) + L(H1, H2, Bv, gv)**2/(4*mu)
        M = max(M, (Phi - A*kap).max())
    return M - (1/3 - b0*(A + 1/3))
def dep_ok(A, Bv, gv, b0):
    # dependent triples: T = -(2B/3) e2(x) - (g m0/3) sum x <= A + C0 + (1-b0) m0^2, x in [-1-m0, 1-m0]^3 (vertices)
    C0 = 1/3 - b0*(A + 1/3); worst = -np.inf
    for m in np.linspace(-1, 1, 81):
        for v in [(1,1,1),(1,1,-1),(1,-1,-1),(-1,-1,-1)]:
            x = np.array(v, float) - m
            e2 = x[0]*x[1] + x[0]*x[2] + x[1]*x[2]
            T = -(2*Bv/3)*e2 - (gv*m/3)*x.sum()
            worst = max(worst, T - (A + C0 + (1 - b0)*m*m))
    return worst
def best_b0(A, Bv, gv):
    lo, hi = 0.0, 0.3
    for _ in range(50):
        mid = (lo + hi)/2
        if worst(A, Bv, gv, mid) <= 0: lo = mid
        else: hi = mid
    return lo
for st in [(3.5, 1.0, 0.0), (3.0, 1.0, -0.5), (3.2, 1.5, -1.0), (2.8, 0.8, -0.3)]:
    r = minimize(lambda p: -best_b0(*p), st, method='Nelder-Mead', options={'xatol': 1e-5, 'fatol': 1e-9, 'maxiter': 3000})
    b = -r.fun; print(np.round(r.x, 4), 'beta0 =', round(b, 6), ' L3 =', round(6/(1 - b), 5), ' dep margin', round(dep_ok(*r.x, b), 4))
print('current: beta0 = 1/13 =', 1/13, ' L3 = 6.5')
