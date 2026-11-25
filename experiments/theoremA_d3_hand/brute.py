# FLOATING POINT brute-force test of (*3): G - B(2X^2 - V) - C X - A kappa <= 1/10 + (9/10) m0^2,
# A = 2, B = 8/9, C = -13 m0/9, triple heights uniform in [c_H - W, c_H + W] or at the chord ends.
import numpy as np
rng = np.random.default_rng(3); N = 2_000_000
m = rng.uniform(-1, 1, N); psi = rng.uniform(0, np.pi/2, N); m0 = rng.uniform(-1, 1, N)
c, s = np.cos(psi), np.sin(psi); cH = m*c; W = np.sqrt(1 - m*m)*s
ends = rng.random((N, 3)) < 0.5; u = np.where(rng.random(N)[:, None] < 0.5, np.where(ends, -1.0, 1.0), rng.uniform(-1, 1, (N, 3)))
t = cH[:, None] + W[:, None]*u
X = t.mean(1) - m0; V = t.var(1); G = (cH - m0)**2 + W*W/2
lhs = G - (8/9)*(2*X*X - V) + (13/9)*m0*X - 2*(1 - c)
print('max LHS - RHS over', N, 'samples:', (lhs - (0.1 + 0.9*m0*m0)).max())
