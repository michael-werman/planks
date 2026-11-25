# FLOATING POINT exploration (round 15): optimise the m0-free scheme of proofs/planar_refinement_hand.md over (A, B, C/m0).
# FLOATING POINT exploration of the m0-free scheme.
# Vertex forms q(c1,c2) = al c1^2 + be c2^2 + 2 ze c1 c2 (c_i = cos theta_i, chord endpoints), delta = theta1-theta2,
# max over sigma: (al+be)/2 + ze cos d + sqrt(((al+be)/2 cos d + ze)^2 + ((al-be)/2)^2 sin^2 d); kappa = 1 - cos(d/2).
import numpy as np
from scipy.optimize import minimize
t = np.linspace(0, 1, 4001)           # t = cos(delta/2)
cd = 2*t*t - 1; sd2 = 4*t*t*(1 - t*t); kap = 1 - t
def qmax(al, be, ze):
    P = (al+be)/2*cd + ze
    return (al+be)/2 + ze*cd + np.sqrt(P*P + ((al-be)/2)**2*sd2)
def forms(B, g, mu):
    # G0 = (3c1^2 + 2c1c2 + 3c2^2)/8 ; split: -B c1 c2 ; L_split = (B - g/2 - 1)(c1 + c2)
    k = B - g/2 - 1
    split = (3/8 + k*k/(4*mu), 3/8 + k*k/(4*mu), 1/8 - B/2 + k*k/(4*mu))
    # (1,1): -B c1^2 ; L_11 = (2B - g - 1) c1 - c2
    l1 = 2*B - g - 1
    one = (3/8 - B + l1*l1/(4*mu), 3/8 + 1/(4*mu), 1/8 - l1/(4*mu))
    return split, one
def worst(A, B, g, b0):
    mu = B - g - b0
    if mu <= 0: return 1e9
    C0 = 1/3 - b0*(A + 1/3)
    s, o = forms(B, g, mu)
    return max((qmax(*s) - A*kap).max(), (qmax(*o) - A*kap).max()) - C0
def best_b0(A, B, g):
    lo, hi = 0.0, 0.2
    for _ in range(60):
        mid = (lo+hi)/2
        if worst(A, B, g, mid) <= 0: lo = mid
        else: hi = mid
    return lo
for st in [(2.94, 1.34, 0.68), (3.0, 1.4, 0.8), (2.9, 1.3, 0.5), (3.2, 1.5, 1.0)]:
    r = minimize(lambda p: -best_b0(*p), st, method='Nelder-Mead', options={'xatol':1e-6,'fatol':1e-10,'maxiter':4000})
    b = -r.fun; print(np.round(r.x, 5), round(b, 6), 'constant', round(2/(1-b), 5))
