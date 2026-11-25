# Blow-up construction (exact integer arithmetic throughout).
# Input: level set L0 subset Z_P (ell elements, every residue a sum of d distinct elements),
#        directions u^(1..ell) in Z^d with every d of them of nonzero determinant.
# Output: S subset Z_N, N = P*prod(a_k), |S| <= ell*T, every residue mod N a sum of d distinct elements.
# usage: python3 blowup.py <name> t [verify]
import sys, itertools, math, time, json, subprocess
from fractions import Fraction
import numpy as np

def det(M):
    M = [list(map(Fraction, r)) for r in M]; n = len(M); D = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None: return Fraction(0)
        if p != c: M[c], M[p] = M[p], M[c]; D = -D
        D *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            for k in range(c, n): M[r][k] -= f * M[c][k]
    return D

def inv(M):
    n = len(M); A = [list(map(Fraction, r)) + [Fraction(int(i == j)) for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0); A[c], A[p] = A[p], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [r[n:] for r in A]

def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]

def lcm(a, b): return a * b // math.gcd(a, b)

def crt(rs, ms):
    x, m = 0, 1
    for r, mi in zip(rs, ms):
        # solve x + m*k = r mod mi
        k = ((r - x) * pow(m, -1, mi)) % mi
        x, m = x + m * k, m * mi
    return x % m

def needed_subsets(L0, P, d, dirs=None):
    """one d-subset of level indices for each residue mod P (first found)."""
    rep = {}
    for J in itertools.combinations(range(len(L0)), d):
        r = sum(L0[j] for j in J) % P
        if dirs is not None and det([[dirs[j][k] for j in J] for k in range(d)]) == 0: continue
        rep.setdefault(r, J)
    assert len(rep) == P, "level set does not cover Z_P"
    return sorted(set(rep.values()))

def build(L0, P, d, dirs, m, c, all_subsets=False, Tover=None):
    ell = len(L0)
    subs = list(itertools.combinations(range(ell), d)) if all_subsets else needed_subsets(L0, P, d, dirs)
    Q = 1
    for J in subs:
        D = det([[dirs[j][k] for j in J] for k in range(d)])
        assert D != 0
        Q = lcm(Q, abs(int(D)))
    Qp = None
    a = [m + ck for ck in c]                           # a_k = m + c_k
    for i, j in itertools.combinations(range(d), 2): assert math.gcd(a[i], a[j]) == 1
    for ak in a: assert math.gcd(ak, P) == 1 and math.gcd(ak, Q) == 1 and (ak - a[0]) % Q == 0
    M = math.prod(a); N = P * M
    # T = max over needed J of max row abs sum of G_J = R^{-1} diag(a) R  (G_J integral)
    T = 0
    for J in subs:
        R = [[dirs[j][k] for j in J] for k in range(d)]      # columns = directions
        G = matmul(matmul(inv(R), [[a[i] if i == j else 0 for j in range(d)] for i in range(d)]), R)
        assert all(x.denominator == 1 for r in G for x in r)
        T = max(T, max(sum(abs(int(x)) for x in r) for r in G))
    Tbound = T
    if Tover is not None: T = Tover
    g = [crt([dirs[j][k] % a[k] for k in range(d)], a) for j in range(ell)]
    X = [sorted(set((i * g[j]) % M for i in range(T))) for j in range(ell)]
    S = []
    level_of = []
    for j in range(ell):
        for x in X[j]:
            S.append(crt([L0[j] % P, x], [P, M])); level_of.append(j)
    return dict(Tbound=Tbound, P=P, d=d, Q=Q, Qp=Qp, a=a, M=M, N=N, T=T, g=g, n=len(S), S=S, level_of=level_of, subs=subs, X=X)

def verify(C, chunk=1 << 22):
    """exact check that the sums over the needed level d-subsets (one element per level, hence distinct
    elements) cover Z_N. Integer arithmetic (int64, all values < 2^62)."""
    N, d = C['N'], C['d']
    assert N * d < 2 ** 62
    hit = np.zeros(N, dtype=np.bool_)
    byl = {}
    for s, j in zip(C['S'], C['level_of']): byl.setdefault(j, []).append(s)
    for J in C['subs']:
        arrs = [np.array(byl[j], dtype=np.int64) for j in J]
        # sums of first d-1 levels
        part = np.zeros(1, dtype=np.int64)
        for A in arrs[:-1]:
            part = ((part[:, None] + A[None, :]) % N).ravel()
        last = arrs[-1]
        step = max(1, chunk // len(last))
        for i in range(0, len(part), step):
            hit[((part[i:i + step, None] + last[None, :]) % N).ravel()] = True
    return bool(hit.all()), int((~hit).sum())

D6 = [(1, 0), (0, 1), (1, 1), (1, -1), (1, 2), (2, 1)]
CONFIGS = {
    # name: (L0, P, d, dirs)
    # name: (L0, P, d, dirs, offsets c_k, admissible m = m0 + step*t)
    'd2_l4': ([0, 1, 2, 4], 6, 2, [(1, 0), (0, 1), (1, 1), (1, -1)], [-1, 1], 6, 6),
    'd2_l12': ([3, 4, 8, 10, 12, 16, 17, 28, 31, 38, 45, 48], 56, 2,
               [D6[c] for c in [0, 1, 2, 0, 2, 1, 0, 3, 4, 0, 4, 5]], [-3, 3], 8, 42),
    'd3_l6b': ([1, 4, 10, 17, 18, 19], 20, 3, [(1, 0, 0), (1, 0, -1), (0, 1, 0), (0, 1, -1), (1, -1, 2), (2, -2, 1)], [-6, 0, 6], 7, 30),
    'd3_l6': ([1, 4, 10, 17, 18, 19], 20, 3, [(0, 0, 1), (1, 2, 2), (1, 1, 1), (1, 0, 2), (1, -1, 1), (0, 2, 1)], [-6, 0, 6], 7, 30),
}

if __name__ == '__main__':
    name, t = sys.argv[1], int(sys.argv[2])
    L0, P, d, dirs, c, m0, step = CONFIGS[name]
    Tover = int(sys.argv[4]) if len(sys.argv) > 4 else None
    t0 = time.time()
    C = build(L0, P, d, dirs, m0 + step * t, c, Tover=Tover)
    out = {k: C[k] for k in ['P', 'd', 'Q', 'a', 'M', 'N', 'Tbound', 'T', 'g', 'n']}
    out['density N/n^d'] = C['N'] / C['n'] ** d
    out['limit P/l^d'] = P / len(L0) ** d
    if len(sys.argv) > 3:
        ok, miss = verify(C)
        out['verified'] = ok; out['missing'] = miss
    out['runtime_s'] = round(time.time() - t0, 2)
    print(json.dumps(out))
