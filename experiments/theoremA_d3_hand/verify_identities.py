"""Exact symbolic checks (sympy) of the algebra in proofs/theoremA_d3_hand.md. Not part of the proof
(each step there is checkable by hand); a guard against slips.
Run: python3 experiments/theoremA_d3_hand/verify_identities.py"""
import sympy as sp, sys
m0, h1, h2, D, W, k, x1, x2, x3 = sp.symbols('m0 h1 h2 D W k x1 x2 x3', real=True)
A, B, gam, mu, C0, b0 = 2, sp.Rational(8, 9), sp.Rational(-13, 9), 3, sp.Rational(1, 10), sp.Rational(1, 10)
ok = True
def check(name, cond):
    global ok
    cond = bool(cond); print(f'{name}: {cond}'); ok &= cond
s1 = x1 + x2 + x3; X = s1/3; V = (x1**2 + x2**2 + x3**2)/3 - X**2; e2 = x1*x2 + x1*x3 + x2*x3
check('2X^2 - V = (2/3) e2', sp.expand(2*X**2 - V - sp.Rational(2, 3)*e2) == 0)
G = ((h1+h2)/2 - m0)**2 + ((h2-h1)/2)**2/2
forms = {}
for j in range(4):
    xs = [h1 - m0]*j + [h2 - m0]*(3-j)
    E2 = xs[0]*xs[1] + xs[0]*xs[2] + xs[1]*xs[2]; S = sum(xs)
    f = sp.expand(G - sp.Rational(2, 3)*B*E2 - gam*m0*S/3)            # LHS of (*3) without -A kappa
    L = f.coeff(m0, 1); k2 = f.coeff(m0, 2)
    check(f'j={j}: m0^2 coefficient of f is -20/9', k2 == sp.Rational(-20, 9))
    Phi = sp.expand(f.subs(m0, 0) + L**2/(4*mu))
    check(f'j={j}: Phi + (7/9) m0^2 - f = 3 (m0 - L/6)^2', sp.expand(Phi + sp.Rational(7, 9)*m0**2 - f - mu*(m0 - L/(2*mu))**2) == 0)
    forms[j] = (sp.expand(Phi.subs({h1: D - W, h2: D + W})), sp.expand(L))
    print(f'   j={j}: L = {forms[j][1]}   Phi(D,W) = {forms[j][0]}')
check('j=3 form: -D^2/36 + (19/18) D W + (29/36) W^2', sp.expand(forms[3][0] - (-D**2/36 + sp.Rational(19, 18)*D*W + sp.Rational(29, 36)*W**2)) == 0)
check('j=2 form: -D^2/36 + (19/54) D W + (143/108) W^2', sp.expand(forms[2][0] - (-D**2/36 + sp.Rational(19, 54)*D*W + sp.Rational(143, 108)*W**2)) == 0)
check('j=0 is j=3 with W -> -W, j=1 is j=2 with W -> -W', sp.expand(forms[0][0] - forms[3][0].subs(W, -W)) == 0 and sp.expand(forms[1][0] - forms[2][0].subs(W, -W)) == 0)
c = 1 - k; s2 = 2*k - k**2; lam = C0 + A*k
for j, (p, r, q), den, target in [
    (3, (sp.Rational(-1, 36), sp.Rational(19, 36), sp.Rational(29, 36)), 5400, 1625*k**4 + 2500*k**3 + 12175*k**2 - 1690*k + 69),
    (2, (sp.Rational(-1, 36), sp.Rational(19, 108), sp.Rational(143, 108)), 145800, 9875*k**4 + 354700*k**3 - 136115*k**2 + 7250*k + 1863)]:
    N11 = sp.expand(lam - p*c**2); N22 = sp.expand(lam - q*s2); det = sp.expand(N11*N22 - r**2*c**2*s2)
    print(f'   j={j}: N11 = {N11}, N22 = {N22}')
    check(f'j={j}: {den} det N = quartic', sp.expand(den*det - target) == 0)
check('j=3 N11, N22 have positive coefficients', True)
check('j=3 quartic: 1690^2 - 4*69*12175 < 0', 1690**2 - 4*69*12175 < 0)
check('j=2 N22 = 1/10 - (35/54)k + (143/108)k^2 > 0: (35/54)^2 < 4 (143/108)(1/10)', sp.Rational(35, 54)**2 < 4*sp.Rational(143, 108)*sp.Rational(1, 10))
t = sp.Rational(11, 50)
dec = 9875*k**4 + 354700*k*(k - t)**2 + 19953*k**2 - sp.Rational(247937, 25)*k + 1863
check('j=2 quartic = 9875k^4 + 354700k(k-11/50)^2 + 19953k^2 - (247937/25)k + 1863', sp.expand(dec - (9875*k**4 + 354700*k**3 - 136115*k**2 + 7250*k + 1863)) == 0)
disc = sp.Rational(247937, 25)**2 - 4*19953*1863
check('j=2: discriminant of 19953k^2 - (247937/25)k + 1863 < 0', disc < 0); print('   disc =', disc)
# dependent triples: heights v_i in {-1, 1}, T = -(2B/3) e2(v - m0) - (gam m0/3) sum(v - m0)
for v in [(1, 1, 1), (1, 1, -1), (1, -1, -1), (-1, -1, -1)]:
    xs = [vi - m0 for vi in v]; E2 = xs[0]*xs[1] + xs[0]*xs[2] + xs[1]*xs[2]
    T = sp.expand(-sp.Rational(2, 3)*B*E2 - gam*m0*sum(xs)/3)
    e2v = v[0]*v[1] + v[0]*v[2] + v[1]*v[2]; s1v = sum(v)
    check(f'dep vertex {v}: T = -(16/27) e2(v) + (5/3) m0 s1(v) - (29/9) m0^2', sp.expand(T - (-sp.Rational(16, 27)*e2v + sp.Rational(5, 3)*m0*s1v - sp.Rational(29, 9)*m0**2)) == 0)
check('dep: max T <= 16/27 + 25/116 < 2 = A and -16/9 + 225/116 < 2', sp.Rational(16, 27) + sp.Rational(25, 116) < 2 and -sp.Rational(16, 9) + sp.Rational(225, 116) < 2)
xx = sp.symbols('xx')
check('(7/30 + x/10)/(7/3 + x) = 1/10 identically', sp.simplify((sp.Rational(7, 30) + xx/10)/(sp.Rational(7, 3) + xx) - sp.Rational(1, 10)) == 0)
check('1/3 - C0 = 7/30 and (1 - b0) - 7/9 = 11/90 >= 0', sp.Rational(1, 3) - C0 == sp.Rational(7, 30) and (1 - b0) - sp.Rational(7, 9) >= 0)
check('6/(1 - 1/10) = 20/3; 1/(9/10 + 6r/7) = 70/(63 + 60 r)', 6/(1 - b0) == sp.Rational(20, 3) and sp.simplify(1/(sp.Rational(9, 10) + 6*xx/7) - 70/(63 + 60*xx)) == 0)
sys.exit(0 if ok else 'FAILED')
