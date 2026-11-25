"""Generate the direction table used by cert.c.
Entry m (0 <= m <= N) is the integer vector (C_m, S_m) = round(2^62 * (cos(m*pi/N), sin(m*pi/N))),
computed with mpmath at 160-bit precision.  Rounding error per component <= 1/2 (+2^-90), so the true
angle of (C_m, S_m) differs from m*pi/N by less than 2^-61 rad (far below one table unit pi/N ~ 3.4e-6 rad).
Output: table_N.bin  (int64 little-endian, 2*(N+1) values)."""
import sys, struct, mpmath, time
from multiprocessing import Pool

N = int(sys.argv[1]) if len(sys.argv) > 1 else 921600
SC = 2 ** 62

def chunk(rng):
    mpmath.mp.prec = 160
    pi = mpmath.pi
    out = []
    for m in range(*rng):
        c, s = mpmath.cos_sin(m * pi / N)
        C = int(mpmath.nint(c * SC)); S = int(mpmath.nint(s * SC))
        out.append((C, S))
    return out

if __name__ == '__main__':
    t0 = time.time()
    step = 20000
    rngs = [(a, min(a + step, N + 1)) for a in range(0, N + 1, step)]
    with Pool(4) as p:
        res = p.map(chunk, rngs)
    with open(f'table_{N}.bin', 'wb') as f:
        for part in res:
            for C, S in part:
                assert -SC <= C <= SC and 0 <= S <= SC
                f.write(struct.pack('<qq', C, S))
    print('N', N, 'entries', N + 1, 'time %.1fs' % (time.time() - t0))
