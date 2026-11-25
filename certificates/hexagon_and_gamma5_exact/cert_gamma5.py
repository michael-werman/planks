"""CERTIFICATE (exact rational arithmetic, Python Fractions/ints only in the checks).
(1) Upper bound: an explicit rational 5-point configuration whose 10 spanned directions, in the
    stated cyclic order, have consecutive gaps with tan(gap) <= t0 (exact).  Since every gap is in
    (0, pi/2) and 10*arctan(t0) < 2*pi, the chain winds exactly once, so the largest gap is
    <= arctan(t0).  Hence gamma_5 <= arctan(t0) and lambda_2(5) <= sin(arctan(t0)/2) <= S_up.
(2) Root isolation: the minimal polynomial p_S (degree 24) of the numerical optimum's sin(gamma/2)
    changes sign on [S_lo, S_hi] (exact evaluation), and S_up - S_hi is tiny.
The configuration's coordinates are rounded (to 1e-80) from a 4000-digit Newton solution."""
import pickle, time, sys
from fractions import Fraction as Fr
from math import isqrt
t_start = time.time()
xs = pickle.load(open('hp5_4000.pkl', 'rb'))
D = 10**80
def rat(s):
    import mpmath as mp; mp.mp.dps = 200
    return Fr(int(mp.nint(mp.mpf(s) * D)), D)
x2, y2, x3, y3, x4, y4 = [rat(s) for s in xs[:6]]
P = [(Fr(0), Fr(0)), (Fr(1), Fr(0)), (x2, y2), (x3, y3), (x4, y4)]
order = [(0, 1), (2, 3), (0, 2), (3, 4), (0, 4), (1, 2), (1, 4), (2, 4), (0, 3), (1, 3)]
assert len(set(order)) == 10
def vec(c):
    i, j = c; return (P[j][0] - P[i][0], P[j][1] - P[i][1])
tans = []
for m in range(10):
    v = vec(order[m]); w = vec(order[(m + 1) % 10])
    cr = v[0] * w[1] - v[1] * w[0]; dt = v[0] * w[0] + v[1] * w[1]
    # directions are lines (mod pi): flip w if needed so that dot>0 is checked on the representative
    if dt < 0: cr, dt = -cr, -dt
    assert dt > 0 and cr > 0, (m, cr, dt)
    tans.append(cr / dt)
tmax = max(tans)
# t0: a short rational >= tmax
t0 = Fr(-((-tmax.numerator * 10**60) // tmax.denominator), 10**60)
assert t0 >= tmax
# 10*arctan(t0) < 2*pi  <=  arctan(t0) < pi/5  <=  t0 < tan(pi/5)=0.7265...: check t0 < 7/10
assert t0 < Fr(7, 10)
# S_up >= sin(arctan(t0)/2) = sqrt((1 - 1/sqrt(1+t0^2))/2).  Exact: c = 1/sqrt(1+t0^2) >= c_lo,
# then S^2 = (1-c)/2 <= (1-c_lo)/2 and S <= sqrt of that rounded up.
E = 10**70
def sqrt_up(q):   # rational >= sqrt(q)
    a = isqrt(q.numerator * E * E // q.denominator) + 1
    return Fr(a, E)
def sqrt_dn(q):
    a = isqrt(q.numerator * E * E // q.denominator)
    return Fr(a, E)
c_lo = 1 / sqrt_up(1 + t0 * t0)
S_up = sqrt_up((1 - c_lo) / 2)
assert S_up * S_up >= (1 - c_lo) / 2
# (2) root isolation for p_S
pS = [16777216, 0, -100663296, 0, 269484032, 0, -424673280, 0, 437649408, 0, -309854208, 0, 154034176,
      0, -53944320, 0, 13123072, 0, -2137088, 0, 215680, 0, -11392, 0, 193]
def ev(c, x):
    r = Fr(0)
    for a in c: r = r * x + a
    return r
S_lo = Fr(1704415526482740379639763544834119968428560732671000730770678952617145, 10**70)
S_hi = S_lo + Fr(1, 10**70)
sl, sh = ev(pS, S_lo), ev(pS, S_hi)
assert sl * sh < 0
print('max tan(gap) of rational configuration:', float(tmax))
print('t0 =', t0)
print('certified: gamma_5 <= arctan(t0); lambda_2(5) <= S_up =', S_up, '~', float(S_up))
print('p_S has a root in [S_lo, S_hi] =', S_lo, '+ [0, 1e-70]; sign values', float(sl), float(sh))
print('S_up - S_hi =', float(S_up - S_hi))
print('runtime %.2fs' % (time.time() - t_start))
