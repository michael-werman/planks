# FLOATING POINT cross-check (round 15) of the scheme against the original per-pair inequality on a grid.
# FLOATING POINT cross-check of the scheme against the original per-pair inequality
# g_Y = d^2 + W^2/2 - B(d-Y)^2 + B(W-|Y|)^2 - C(d-Y) - A kappa <= C0(m0) for Y in {-W,0,W}, C = g*m0,
# C0(m0) = 1/3 - b0(A+1/3) + (1-b0) m0^2  (which gives beta = b0 exactly at I = 1/3 + m0^2).
import numpy as np
A, B, g, b0 = 2.90812, 1.32757, -0.72097, 0.0792
phi = np.linspace(0, np.pi, 1201)[:, None]; psi = np.linspace(0, np.pi/2, 1201)[None, :]
c, s = np.cos(psi), np.sin(psi); kap = 1 - c
worst = -1
for m0 in np.linspace(-1, 1, 81):
    d = np.cos(phi)*c - m0; W = np.sin(phi)*s; C = g*m0
    v = max((d*d + W*W/2 - B*(d-Y)**2 + B*(W-np.abs(Y))**2 - C*(d-Y) - A*kap).max() for Y in (-W, 0*W, W))
    C0 = 1/3 - b0*(A + 1/3) + (1 - b0)*m0*m0
    worst = max(worst, v - C0)
print('max over m0 grid of (max g_Y - C0(m0)):', worst)
