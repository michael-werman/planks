# FLOATING POINT (forms computed exactly by sympy): d = 3 vertex forms Phi = p D^2 + 2 r D W + q W^2
import numpy as np, sympy as sp
from fractions import Fraction as Fr
h1, h2, m0, B, g, mu, D, W = sp.symbols('h1 h2 m0 B g mu D W')
G = ((h1+h2)/2 - m0)**2 + ((h2-h1)/2)**2/2
def forms(Bv, gv, muv):
    out = []
    for j in (3, 2):
        x1, x2 = h1 - m0, h2 - m0
        e2 = sp.binomial(j, 2)*x1**2 + j*(3-j)*x1*x2 + sp.binomial(3-j, 2)*x2**2
        f = sp.expand((G - sp.Rational(2, 3)*B*e2 - g*m0*(j*x1 + (3-j)*x2)/3).subs({B: Bv, g: gv}))
        Phi = sp.expand((f.subs(m0, 0) + f.coeff(m0, 1)**2/(4*muv)).subs({h1: D - W, h2: D + W}))
        out.append((Phi.coeff(D, 2), Phi.coeff(D, 1).coeff(W, 1)/2, Phi.coeff(W, 2), sp.expand(f.coeff(m0, 2))))
    return out
k = np.linspace(0, 1, 200001); c = 1 - k; s2 = 1 - c*c; cs = c*np.sqrt(s2)
def C0needed(A, fs):
    M = -np.inf
    for p, r, q, _ in fs:
        p, r, q = float(p), float(r), float(q)
        a, d, b = p*c*c, q*s2, r*cs
        lam = (a + d)/2 + np.sqrt(((a - d)/2)**2 + b*b)
        M = max(M, (lam - A*k).max())
    return M
def beta(A, Bv, gv, muv):
    fs = forms(Bv, gv, muv)
    C0 = C0needed(A, fs)
    b0 = (Fr(1, 3) - Fr(C0).limit_denominator(10**12))/(A + Fr(1, 3))
    ok = all(1 - float(b0) - float(k2) >= muv for *_, k2 in fs)
    return float(b0), C0, ok, fs
print('optimum-ish:', beta(2.0115, 0.9017, -1.4071, 2.0*0.9017 + 1.4071 - 0.103)[:3])
for A in [Fr(2), Fr(201, 100)]:
    for Bv in [Fr(9, 10), Fr(8, 9), Fr(11, 12)]:
        for gv in [Fr(-7, 5), Fr(-4, 3), Fr(-3, 2), Fr(-13, 9)]:
            for muv in [Fr(3), Fr(31, 10)]:
                b0, C0, ok, fs = beta(A, sp.Rational(Bv.numerator, Bv.denominator), sp.Rational(gv.numerator, gv.denominator), sp.Rational(muv.numerator, muv.denominator))
                if ok: print(A, Bv, gv, muv, ' beta0 = %.6f  L3 = %.5f  C0 = %.6f' % (b0, 6/(1-b0), C0), [(f[0], f[1], f[2]) for f in fs])
