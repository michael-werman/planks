"""Exact symbolic checks (sympy, rational arithmetic) of the algebra in proofs/planar_refinement_hand.md.
Not part of the proof (every step there is checkable by hand); a guard against slips.
Run: python3 experiments/planar_hand/verify_identities.py"""
import sympy as sp, sys
m0, h1, h2, a, b, D, W, c, s, k, t, phi, psi = sp.symbols('m0 h1 h2 a b D W c s k t phi psi', real=True)
A, B, gam, mu = sp.Rational(35, 12), sp.Rational(4, 3), sp.Rational(-2, 3), sp.Rational(16, 9)
C0 = sp.Rational(27, 352)
ok = True
def check(name, cond):
    global ok
    cond = bool(cond); print(f'{name}: {cond}'); ok &= cond
# chord: endpoint heights h1, h2; G = d^2 + W^2/2 with d = (h1+h2)/2 - m0, W = (h2-h1)/2
G = ((h1+h2)/2 - m0)**2 + ((h2-h1)/2)**2/2
# pair at heights a, b on the chord: X = (a+b)/2 - m0, Z = (a-b)/2; f = G - B X^2 + B Z^2 - C X  (kappa term separate)
X = (a+b)/2 - m0; Z = (a-b)/2; C = gam*m0
f = G - B*X**2 + B*Z**2 - C*X
check('B X^2 - B Z^2 = B x_a x_b', sp.expand(B*X**2 - B*Z**2 - B*(a-m0)*(b-m0)) == 0)
# completing the square in m0: f = f(0) + m0 L + (1 - B + gam) m0^2, L = (B - gam/2)(a+b) - (h1+h2)
L = (B - gam/2)*(a+b) - (h1+h2)
check('f = f|_{m0=0} + m0 L - m0^2', sp.expand(f - (f.subs(m0, 0) + m0*L - m0**2)) == 0)
Phi = f.subs(m0, 0) + L**2/(4*mu)
check('Phi + (7/9) m0^2 - f = mu (m0 - L/(2 mu))^2 >= 0', sp.expand(Phi + sp.Rational(7, 9)*m0**2 - f - mu*(m0 - L/(2*mu))**2) == 0)
# split vertex (a, b) = (h1, h2)
Ps = sp.expand(Phi.subs({a: h1, b: h2}))
check('split: Phi = -(h1+h2)^2/48 + (11/24)(h1-h2)^2', sp.expand(Ps - (-(h1+h2)**2/48 + sp.Rational(11, 24)*(h1-h2)**2)) == 0)
# same-endpoint vertex (a, b) = (h1, h1), with h1 = D - W, h2 = D + W
P1 = sp.expand(Phi.subs({a: h1, b: h1}).subs({h1: D - W, h2: D + W}))
check('same end: Phi = -D^2/12 + (17/12) D W + (35/48) W^2', sp.expand(P1 - (-D**2/12 + sp.Rational(17, 12)*D*W + sp.Rational(35, 48)*W**2)) == 0)
P2 = sp.expand(Phi.subs({a: h2, b: h2}).subs({h1: D - W, h2: D + W}))
check('other end: same with W -> -W', sp.expand(P2 - P1.subs(W, -W)) == 0)
# split bound: (11/6)(1 - c^2) - A(1 - c) <= C0 for all c, i.e. C0 - that = (11/6)(k - 9/44)^2 with k = 1 - c
check('split: C0 - [(11/6)(1-c^2) - A(1-c)] = (11/6)(1 - c - 9/44)^2', sp.expand(C0 - (sp.Rational(11, 6)*(1-c**2) - A*(1-c)) - sp.Rational(11, 6)*(1 - c - sp.Rational(9, 44))**2) == 0)
# eigenvalue: matrix of -D^2/12 + (17/12) D W + (35/48) W^2 in (cos phi, sin phi), D = c cos phi, W = s sin phi
M = sp.Matrix([[-c**2/12, sp.Rational(17, 24)*c*s], [sp.Rational(17, 24)*c*s, sp.Rational(35, 48)*s**2]])
lam = (35*s**2 - 4*c**2 + sp.sqrt((35*s**2 + 4*c**2)**2 + 4624*c**2*s**2))/96
check('lambda_max is a root of det(M - lambda I)', sp.simplify((M - lam*sp.eye(2)).det()) == 0)
# quartic: with s^2 = 1 - c^2, c = 1 - k, R = 96(C0 + A k) - 35 s^2 + 4 c^2, S = (35 s^2 + 4 c^2)^2 + 4624 c^2 s^2
cc = 1 - k; ss2 = 1 - cc**2
R = 96*(C0 + A*k) - 35*ss2 + 4*cc**2
S = (35*ss2 + 4*cc**2)**2 + 4624*cc**2*ss2
target = 627264*k**4 + 133584*k**3 + 7406938*k**2 - 623524*k + 13689
check('121 (R^2 - S) = Q(k)', sp.expand(121*(R**2 - S) - target) == 0)
check('R = 39 k^2 + 202 k + 125/11 (so R > 0 for k >= 0)', sp.expand(R - (39*k**2 + 202*k + sp.Rational(125, 11))) == 0)
disc = 623524**2 - 4*7406938*13689
check('discriminant of 7406938 k^2 - 623524 k + 13689 < 0', disc < 0)
print('discriminant =', disc)
# final constant
I0 = sp.Rational(1, 3); beta0 = (I0 - C0)/(A + I0)
check('beta0 = 271/3432', beta0 == sp.Rational(271, 3432))
check('2/(1 - beta0) = 6864/3161 > 76/35 > 2.1714', 2/(1 - beta0) == sp.Rational(6864, 3161) and sp.Rational(6864, 3161) > sp.Rational(76, 35) > sp.Rational(21714, 10000))
check('beta(m0) = (271/1056 + (2/9) m0^2)/(13/4 + m0^2) is increasing: (2/9)(13/4) > 271/1056', sp.Rational(2, 9)*sp.Rational(13, 4) > sp.Rational(271, 1056))
sys.exit(0 if ok else 'FAILED')
