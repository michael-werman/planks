"""Floating-point exploration (NOT a certificate): j-flats through j+1 points in B^d.

1. Lower-bound density f_a(x) = (1-|x|^2)^(-a): checks numerically that the integral over a j-flat at
   distance h is R^(j-2a) * (|S^{j-1}|/2) B(j/2, 1-a), R = sqrt(1-h^2), and the constant
   c_a(d,j) = Gamma(j/2+1-a) Gamma((d-j)/2+1) / Gamma(d/2+1-a), with limit (a -> 1, j >= 2)
   c(d,j) = Gamma(j/2) Gamma((d-j)/2+1) / Gamma(d/2).
2. Construction: R^d = R^{j-1} x R^m, m = d-j+1. Groups r = 0..j at moment-curve nodes mu(r) in
   R^{j-1}; group r carries vectors b/c_r, b in B_r, where B_0 = {(g, h)}, B_r = {(M^r g, 0)}
   (g in {0..M-1}^{m-1}), shifted to centre the grid; c_r the affine dependency of mu(0..j).
   Checks (a) that the actual flat through a transversal (j+1)-set at small s is close to
   R^{j-1} + R w_F with w_F = sum of the b's, and (b) estimates the covering radius delta of
   {w_F} on S^{m-1} (one block, chart region only) against the bound sqrt(m-1)/(N-1), N = M^(j+1).
Usage: python3 flats_middle_dims.py
"""
import numpy as np, itertools
from math import gamma, pi, sqrt, comb
rng = np.random.default_rng(3)
def S(k): return 2*pi**(k/2)/gamma(k/2)
def B(x, y): return gamma(x)*gamma(y)/gamma(x+y)

print("1. flat integrals of (1-|x|^2)^(-a)")
for j, a in [(2, 0.9), (3, 0.6), (2, 0.5)]:
    for h in [0.0, 0.5, 0.9]:
        R = sqrt(1-h*h); N = 2_000_000
        y = rng.uniform(-R, R, (N, j)); r2 = (y**2).sum(1); ok = r2 < R*R
        est = (2*R)**j * np.mean(np.where(ok, (R*R - np.where(ok, r2, 0) + 1e-300)**(-a), 0))
        exact = R**(j-2*a) * S(j)/2 * B(j/2, 1-a)
        print(f"  j={j} a={a} h={h}: MC {est:.4f}  formula {exact:.4f}")
print("   c(d,j) = Gamma(j/2)Gamma((d-j)/2+1)/Gamma(d/2):")
for d in range(3, 8):
    row = []
    for j in range(1, d):
        c = 1.0 if j == 1 else gamma(j/2)*gamma((d-j)/2+1)/gamma(d/2)
        vol = gamma(j/2+1)*gamma((d-j)/2+1)/gamma(d/2+1)      # plain volume bound
        row.append(f"j={j}: {c:.3f} (vol {vol:.3f})")
    print(f"  d={d}: " + "; ".join(row))

def construction(d, j, M):
    m = d - j + 1
    nodes = np.arange(j+1, dtype=float)
    c = np.array([1/np.prod([r-k for k in range(j+1) if k != r]) for r in range(j+1)])
    mu = np.stack([nodes**p for p in range(1, j)], 1) if j > 1 else np.zeros((j+1, 0))
    N = M**(j+1); h = (N-1)/2
    G = np.array(list(itertools.product(range(M), repeat=m-1)), dtype=float)
    groups = []
    for r in range(j+1):
        if r == 0: Bv = np.hstack([G - (N-1)/2, np.full((len(G), 1), h)])
        else:      Bv = np.hstack([M**r * G, np.zeros((len(G), 1))])
        groups.append(Bv / c[r])
    return m, mu, c, groups, N

print("\n2. construction checks")
for d, j, M in [(4, 2, 3), (5, 2, 2), (5, 3, 2), (4, 1, 4)]:
    m, mu, c, groups, N = construction(d, j, M)
    # (a) actual flat vs predicted limit, small s
    s = 1e-6; kappa = 0.1; errs = []
    for _ in range(50):
        idx = [rng.integers(len(g)) for g in groups]
        Q = np.array([np.concatenate([kappa*mu[r], s*groups[r][idx[r]]]) for r in range(j+1)])
        D = (Q[1:] - Q[0]).T                     # d x j direction vectors
        Qb, _ = np.linalg.qr(D)
        wF = sum(c[r]*groups[r][idx[r]] for r in range(j+1))
        pred = np.zeros((d, j)); 
        for k in range(j-1): pred[k, k] = 1
        pred[j-1:, j-1] = wF/np.linalg.norm(wF)
        # principal angle distance between subspaces
        sv = np.linalg.svd(Qb.T @ pred, compute_uv=False)
        errs.append(np.sqrt(max(0, 1 - sv.min()**2)))
    # (b) covering radius of directions w_F (one block covers the chart |x_k| <= |x_m|)
    W = [np.zeros(m)]
    for r in range(j+1):
        W = [w + c[r]*g for w in W for g in groups[r]]
    W = np.array(W); W /= np.linalg.norm(W, axis=1)[:, None]
    X = rng.normal(size=(200_000, m)); X /= np.linalg.norm(X, axis=1)[:, None]
    X = X[np.all(np.abs(X[:, :-1]) <= np.abs(X[:, -1:]), axis=1)]   # chart of the last coordinate
    best = np.zeros(len(X))
    for k in range(0, len(W), 2000): best = np.maximum(best, np.abs(X @ W[k:k+2000].T).max(1))
    delta = np.arccos(best.min())
    n = m*(j+1)*M**(m-1)
    print(f"  d={d} j={j} M={M}: flat error at s=1e-6 <= {max(errs):.1e}; "
          f"delta ~ {delta:.4f} <= bound {sqrt(m-1)/(N-1):.4f}; n = {n}, "
          f"n^((j+1)/(d-j)) sin(delta) ~ {n**((j+1)/(d-j))*np.sin(delta):.1f}")
