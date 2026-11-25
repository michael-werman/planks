"""Soundness self-test input (FLOATING POINT generation; the check itself is in cert.c trace mode).
Generates configurations with max gap < theta_deg (near stored local minima, plus perturbations),
normalises them exactly as in the certificate (longest pair -> (0,0),(1,0); x-reflection so that
x_min + x_max <= 1; sort by x; y-reflection so y2 >= 0) and prints grid coordinates (units 2^-K)."""
import sys, glob, math, itertools, numpy as np
sys.path.insert(0, '.')
from gamma import maxgap

n = int(sys.argv[1]); theta = math.radians(float(sys.argv[2])); K = int(sys.argv[3]); want = int(sys.argv[4])
rng = np.random.default_rng(int(sys.argv[5]) if len(sys.argv) > 5 else 0)

def normalise(P):
    P = np.array(P, float)
    d = [(np.hypot(*(P[i] - P[j])), i, j) for i, j in itertools.combinations(range(len(P)), 2)]
    _, i, j = max(d)
    v = P[j] - P[i]; s = np.hypot(*v); c, sn = v / s
    R = np.array([[c, sn], [-sn, c]])
    Q = (P - P[i]) @ R.T / s
    others = [k for k in range(len(P)) if k not in (i, j)]
    O = Q[others]
    if O[:, 0].min() + O[:, 0].max() > 1: O[:, 0] = 1 - O[:, 0]
    O = O[np.argsort(O[:, 0])]
    if O[0, 1] < 0: O[:, 1] = -O[:, 1]
    return O

seeds = []
for f in glob.glob(f'gam_n{n}_*.npz') + glob.glob(f'gam_n{n}_s*.npz'):
    for P in np.load(f)['P']:
        if maxgap(P) < theta: seeds.append(P)
out = []
tries = 0
while len(out) < want and tries < 200 * want:
    tries += 1
    P = seeds[rng.integers(len(seeds))]
    scale = 10 ** rng.uniform(-6, -1.5)
    Q = P + rng.normal(size=P.shape) * scale * np.abs(P).max()
    g = maxgap(Q)
    if g < theta - 1e-9:
        out.append(Q)
for Q in out:
    O = normalise(Q)
    G = [int(math.floor(v * 2 ** K)) for v in O.ravel()]
    print(' '.join(map(str, G)))
print(f'# seeds {len(seeds)} generated {len(out)} tries {tries}', file=sys.stderr)
