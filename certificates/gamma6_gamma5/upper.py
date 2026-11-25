"""Rigorous UPPER bound gamma_n <= G(P) for an explicit configuration P with rational coordinates.
Directions via mpmath interval arithmetic (iv.atan2) at 200 bits; the configuration is the stored
numerical optimiser (decimal strings, read as exact rationals).  Prints an upper bound for max gap (deg)."""
import json, sys, itertools, time
from fractions import Fraction
import mpmath
from mpmath import iv, mpf
iv.prec = 200

def upper(points):
    P = [(Fraction(x), Fraction(y)) for x, y in points]
    n = len(P)
    dirs = []
    for i, j in itertools.combinations(range(n), 2):
        dx = P[j][0] - P[i][0]; dy = P[j][1] - P[i][1]
        if dy < 0 or (dy == 0 and dx < 0): dx, dy = -dx, -dy      # direction mod pi in [0, pi)
        a = iv.atan2(iv.mpf(dy.numerator) / dy.denominator, iv.mpf(dx.numerator) / dx.denominator)
        dirs.append(a)
    # sort by midpoint; check intervals are disjoint so the order is certain
    dirs.sort(key=lambda a: a.mid)
    for a, b in zip(dirs, dirs[1:]):
        assert a.b < b.a, 'order not certified'
    gaps = [dirs[k + 1] - dirs[k] for k in range(len(dirs) - 1)] + [dirs[0] + iv.pi - dirs[-1]]
    gmax_hi = max(g.b for g in gaps)
    return gmax_hi * 180 / iv.pi, gaps

if __name__ == '__main__':
    t0 = time.time()
    cfg = json.load(open(sys.argv[1]))
    for key in ('gamma_5_optimum', 'gamma_6_optimum'):
        ub, gaps = upper(cfg[key]['points'])
        print(key, ': max gap <= %s deg (upper end of interval)' % mpmath.nstr(ub.b, 20))
        print('   gaps (deg, sorted):', sorted(round(float((g * 180 / iv.pi).mid), 6) for g in gaps))
    print('runtime %.2fs' % (time.time() - t0))
