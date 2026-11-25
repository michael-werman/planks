"""Full rigorous check of the direction table: for every m, |angle(C_m,S_m) - m*pi/N| < 2^-60 (mpmath iv, 100 bits)."""
import sys, numpy as np, time
from mpmath import iv
iv.prec = 100
N = 921600
T = np.fromfile(f'table_{N}.bin', dtype='<i8').reshape(-1, 2)
assert T.shape == (N + 1, 2)
r, P = int(sys.argv[1]), int(sys.argv[2])
t0 = time.time(); worst = 0; bound = iv.mpf(2) ** -60
for m in range(r, N + 1, P):
    C, S = int(T[m, 0]), int(T[m, 1])
    d = abs(iv.atan2(iv.mpf(S), iv.mpf(C)) - iv.mpf(m) * iv.pi / N)
    assert d.b < bound.a, m
    if d.b > worst: worst = d.b
print(f'part {r}/{P}: all entries OK, max error bound {float(worst):.3e} rad, {time.time()-t0:.0f}s')
