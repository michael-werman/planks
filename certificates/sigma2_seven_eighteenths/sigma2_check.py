"""Exact checks for proofs/sigma2_seven_eighteenths.md (integers and Fractions only).

(1) The seed L ⊂ Z_56 (12 elements) has every residue as a sum of two distinct elements.
(2) Q = lcm(1..11) = 27720 and 56 | Q.
(3) For every pair 1 <= s < t <= 12 the matrix G(s,t) of Lemma 2 has integer entries, its columns
    lie in the kernel, det G = m(m+Q), and its row sums of |entries| are <= m + 23Q.  Entries are
    affine in m and det is quadratic in m, so the identities are checked as polynomials in m
    (exact coefficients).  The row-sum bound is checked for all m >= 1 by the same coefficients.
(5) Lemma 2 at the real parameters Q = 27720, m = 6Q+1 and m = 1000003, all 66 pairs, by the
    arc test (see section (5) below).
(4) Brute-force sanity checks of Lemma 2 and of the full lift for small analogues:
    ell' indices, Q' = lcm(1..ell'-1) (times the seed modulus), several m coprime to Q'.
"""
import sys, time
from fractions import Fraction as F
from math import gcd
from itertools import combinations

def lcm(a, b): return a * b // gcd(a, b)

def restricted_cover(L, P):
    """Return dict r -> (i, j) with i < j indices, L[i]+L[j] = r mod P, or None if some r missed."""
    rep = {}
    for i, j in combinations(range(len(L)), 2):
        rep.setdefault((L[i] + L[j]) % P, (i, j))
    return rep if len(rep) == P else None

t0 = time.time()
ok = True
# (1) seed
L = [0, 1, 12, 15, 22, 29, 32, 43, 44, 48, 50, 52]
P = 56
assert len(set(x % P for x in L)) == 12
rep = restricted_cover(L, P)
print("(1) seed covers Z_56 with distinct summands:", rep is not None)
ok &= rep is not None
if rep:
    print("    representations r: (l_s, l_t):",
          " ".join(f"{r}:({L[i]},{L[j]})" for r, (i, j) in sorted(rep.items())))

# (2) Q
Q = 1
for k in range(1, 12): Q = lcm(Q, k)
print("(2) Q = lcm(1..11) =", Q, "; 56 | Q:", Q % 56 == 0)
ok &= (Q == 27720 and Q % 56 == 0)

# (3) G(s,t) as polynomials in m: entry = c1*m + c0
def G_poly(s, t, Q):
    d = t - s
    return [[(1, F(-Q * s, d)), (0, F(-Q * t, d))],
            [(0, F(Q * s, d)), (1, F(Q * t, d))]]

ok3 = True
for s, t in combinations(range(1, 13), 2):
    G = G_poly(s, t, Q)
    for row in G:
        for c1, c0 in row:
            ok3 &= F(c0).denominator == 1
    # columns in kernel: x+y = m*(something) exactly m ; s x + t y = multiple of (m+Q)
    for col in range(2):
        x1, x0 = G[0][col]; y1, y0 = G[1][col]
        # x + y as polynomial must be exactly m (≡ 0 mod a = m)
        ok3 &= (x1 + y1, x0 + y0) == (1, 0)
        # s x + t y must be (s or t)*(m + Q): coefficient of m equals constant/Q
        a1, a0 = s * x1 + t * y1, s * x0 + t * y0
        ok3 &= a0 == a1 * Q and F(a1).denominator == 1
    # det = (c1 m + c0)(...) - (...)(...) as quadratic: compare coefficients with m^2 + Q m
    (p1, p0), (q1, q0) = G[0]; (r1, r0), (u1, u0) = G[1]
    det = (p1 * u1 - q1 * r1, p1 * u0 + p0 * u1 - q1 * r0 - q0 * r1, p0 * u0 - q0 * r0)
    ok3 &= det == (1, Q, 0)
    # row sums of |entries| <= m + 23Q for all m >= 1: |c1 m + c0| <= m + |c0| (c1 in {0,1})
    for row in G:
        bound_const = sum(abs(c0) for c1, c0 in row)
        ok3 &= bound_const <= 23 * Q and sum(c1 for c1, c0 in row) == 1
print("(3) G(s,t): integer, kernel, det = m(m+Q), row sums <= m + 23Q, all 66 pairs:", ok3)
ok &= ok3

# (4a) Lemma 2 brute force, small analogue
def lemma2_brute(ellp, Qp, m, Tbound):
    a, b = m, m + Qp
    M = a * b
    # CRT element g_t: ≡ 1 mod a, ≡ t mod b
    inv_a_mod_b = pow(a, -1, b)
    def crt(x, y):  # z ≡ x mod a, ≡ y mod b
        return (x + a * ((y - x) * inv_a_mod_b % b)) % M
    res = True; worst = 0
    for s, t in combinations(range(1, ellp + 1), 2):
        gs, gt = crt(1, s), crt(1, t)
        Xs = [(i * gs) % M for i in range(Tbound)]
        Xt = [(i * gt) % M for i in range(Tbound)]
        cov = bytearray(M)
        for x in Xs:
            for y in Xt:
                cov[(x + y) % M] = 1
        res &= all(cov)
    return res

ok4 = True
for ellp in (4, 5):
    Qp = 1
    for k in range(1, ellp): Qp = lcm(Qp, k)
    for m in [x for x in range(5, 60) if gcd(x, Qp) == 1][:3] + [x for x in (101, 203, 301, 401) if gcd(x, Qp) == 1]:
        T = m + (2 * ellp - 1) * Qp
        r = lemma2_brute(ellp, Qp, m, T)
        ok4 &= r
        print(f"(4a) ell'={ellp} Q'={Qp} m={m} M={m*(m+Qp)} T={T} (nontrivial: T < M: {T < m*(m+Qp)}): X_s+X_t = Z_M for all pairs: {r}")

# (4b) full lift, toy seed: a 5-subset of Z_9 (p_2(5) = 9)
def find_seed(k, P):
    for rest in combinations(range(1, P), k - 1):
        Ls = [0, *rest]
        r = restricted_cover(Ls, P)
        if r: return Ls, r
    return None, None

Ls, reps = find_seed(5, 9)
print("(4b) toy seed", Ls, "covers Z_9:", Ls is not None)
ok4 &= Ls is not None
ellp, Pp = 5, 9
Qp = lcm(lcm(lcm(1, 2), lcm(3, 4)), Pp)  # 36
for m in [x for x in range(5, 60) if gcd(x, Qp) == 1][:2] + [301, 401]:
    a, b = m, m + Qp; M = a * b; N = Pp * M
    T = m + (2 * ellp - 1) * Qp
    # element of Z_N from (r mod P, z mod M): CRT
    inv = pow(Pp, -1, M)
    def lift(r, z): return (r + Pp * ((z - r) * inv % M)) % N
    inv_a_mod_b = pow(a, -1, b)
    def crt(x, y): return (x + a * ((y - x) * inv_a_mod_b % b)) % M
    S = set()
    for t in range(1, ellp + 1):
        gt = crt(1, t)
        for i in range(T):
            S.add(lift(Ls[t - 1], (i * gt) % M))
    S = sorted(S)
    cov = bytearray(N)
    for i in range(len(S)):
        si = S[i]
        for j in range(i + 1, len(S)):
            cov[(si + S[j]) % N] = 1
    r = all(cov)
    ok4 &= r
    print(f"(4b) lift m={m} (nontrivial: T < M: {T < M}): N={N}, |S|={len(S)} <= {ellp*T}: S+S (distinct) = Z_N: {r}")

ok &= ok4

# (5) Lemma 2 at the real parameters, by the arc test (exact integers).
# X_s + X_t = Z_M  <=>  g_t^{-1}(X_s + X_t) = {i*u + j : 0 <= i, j < T} = Z_M, u = g_s * g_t^{-1},
# <=> the arcs [i*u, i*u + T) (mod M), 0 <= i < T, cover Z_M
# <=> the sorted start points i*u mod M have all circular gaps <= T.
ok5 = True
for m in (6 * Q + 1, 1000003):
    assert gcd(m, Q) == 1
    a, b = m, m + Q; M = a * b; T = m + 23 * Q
    inv_a_mod_b = pow(a, -1, b)
    def crt(x, y): return (x + a * ((y - x) * inv_a_mod_b % b)) % M
    g = {t: crt(1, t) for t in range(1, 13)}
    worst = 0
    for s_, t_ in combinations(range(1, 13), 2):
        u = g[s_] * pow(g[t_], -1, M) % M
        st = sorted(i * u % M for i in range(T))
        gaps = max(max(st[i + 1] - st[i] for i in range(T - 1)), st[0] + M - st[-1])
        worst = max(worst, gaps)
        ok5 &= gaps <= T
    print(f"(5) real Q: m={m} M={M} T={T}: all 66 pairs cover Z_M (max circular gap {worst} <= T): {worst <= T}")
ok &= ok5
dt = time.time() - t0
print(f"ALL CHECKS PASSED: {ok}   runtime {dt:.1f} s")
sys.exit(0 if ok else 1)
