"""Independent re-check of sampled pruned leaves of cert.c.
Written separately from cert.c: direction intervals from mpmath interval arithmetic (iv.atan2 at 100 bits,
no table), converted to integer units (NT units = pi) by outward floor/ceil; domain test with Python ints;
tests P, J, M re-implemented.  Narrower intervals can only make P/J/M prune more, so every leaf pruned
by cert.c must also be pruned here.  usage: recheck.py n THETA_units K leaves_file [max]"""
import sys, itertools, math, time
from mpmath import iv
iv.prec = 100
NT = 921600

def main():
    n = int(sys.argv[1]); TH = int(sys.argv[2]); K = int(sys.argv[3]); fn = sys.argv[4]
    mx = int(sys.argv[5]) if len(sys.argv) > 5 else 10 ** 9
    ONE = 1 << K; pairs = list(itertools.combinations(range(n), 2)); L = len(pairs); E = L * TH - NT
    PIU = iv.pi

    def units_interval(lo_rad, hi_rad):
        a = iv.mpf(lo_rad) * NT / PIU; b = iv.mpf(hi_rad) * NT / PIU
        return math.floor(a.a), math.ceil(b.b)

    def ang(x, y):          # interval enclosure of atan2(y, x) (radians), exact integer inputs
        return iv.atan2(iv.mpf(y), iv.mpf(x))

    def direction(ax0, ax1, ay0, ay1):
        """(lo, w) in units with all directions mod pi in [lo, lo+w], or None if unknown"""
        if ax0 <= 0 <= ax1 and ay0 <= 0 <= ay1: return None
        corners = [(ax0, ay0), (ax1, ay0), (ax0, ay1), (ax1, ay1)]
        # flip into the half plane y>0 or x>0 (direction mod pi unchanged); angles then in (-pi/2, pi)
        if ay1 < 0 or (ay0 <= 0 <= ay1 and ax1 < 0):
            corners = [(-x, -y) for x, y in corners]
        A = [ang(x, y) for x, y in corners]
        lo = min(a.a for a in A); hi = max(a.b for a in A)
        ul, uh = units_interval(lo, hi)
        if uh - ul >= NT: return None
        return ul % NT, uh - ul

    def domain_bad(pts):
        (x0, x1, y0, y1) = zip(*pts)
        if pts[2][3] < 0: return True
        for k in range(2, n - 1):
            if pts[k][0] > pts[k + 1][1]: return True
        if pts[2][0] + pts[n - 1][0] > ONE: return True
        for i, j in pairs:
            a0 = pts[j][0] - pts[i][1]; a1 = pts[j][1] - pts[i][0]; b0 = pts[j][2] - pts[i][3]; b1 = pts[j][3] - pts[i][2]
            mxx = a0 if a0 > 0 else (-a1 if a1 < 0 else 0); myy = b0 if b0 > 0 else (-b1 if b1 < 0 else 0)
            if mxx * mxx + myy * myy > ONE * ONE: return True
        return False

    def dirs_prune(pts):
        D = []
        for i, j in pairs:
            D.append(direction(pts[j][0] - pts[i][1], pts[j][1] - pts[i][0], pts[j][2] - pts[i][3], pts[j][3] - pts[i][2]))
        unk = sum(d is None for d in D)
        for p in range(L):
            if D[p] is None: continue
            lp, wp = D[p]
            R = {}
            for q in range(L):
                if q == p: continue
                if D[q] is None: R[q] = (1, L - 1); continue
                lq, wq = D[q]; s = (lq - lp) % NT
                if not (wp < s and s + wq < NT):        # overlapping arcs
                    h = min(max(wp, s + wq), max(wq, NT - s + wp))
                    if h <= TH - E: return 'P'
                    R[q] = (1, L - 1); continue
                lmin, lmax = s - wp, s + wq
                jmin = 0; jmax = unk
                for r in range(L):
                    if r in (p, q) or D[r] is None: continue
                    t = (D[r][0] - lp) % NT; wr = D[r][1]
                    if t > wp and t + wr < s: jmin += 1
                    if t < s + wq or t + wr > NT: jmax += 1
                if not any(lmin < m * TH and lmax > m * TH - E for m in range(jmin + 1, jmax + 2)): return 'J'
                a = max(1, lmin // TH + 1); b = min(L - 1, -((-(lmax + E)) // TH) - 1)
                if a > b: return 'J'
                R[q] = (a, b)
            # perfect matching lines <-> ranks 1..L-1 with interval rank sets: Hall's condition on intervals
            for lo_ in range(1, L):
                for hi_ in range(lo_, L):
                    inside = sum(1 for (a, b) in R.values() if lo_ <= a and b <= hi_)
                    if inside > hi_ - lo_ + 1: return 'M'
                    cover = sum(1 for (a, b) in R.values() if not (b < lo_ or a > hi_))
                    if cover < hi_ - lo_ + 1: return 'M'
        return None

    t0 = time.time(); cnt = 0; bad = 0; codes = {}
    for line in open(fn):
        v = list(map(int, line.split())); code, X = v[0], v[1:]
        pts = [(0, 0, 0, 0), (ONE, ONE, 0, 0)] + [tuple(X[4 * k:4 * k + 4]) for k in range(n - 2)]
        res = 'D' if domain_bad(pts) else dirs_prune(pts)
        cnt += 1; codes[res] = codes.get(res, 0) + 1
        if res is None:
            bad += 1; print('NOT PRUNED by recheck:', line.strip())
        if cnt >= mx: break
    print(f'rechecked {cnt} leaves from {fn}: python verdicts {codes}; failures {bad}; {time.time()-t0:.1f}s')

main()
