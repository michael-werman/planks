# Finite seeds (found by heuristic search; verified here exactly with Python ints by enumerating all
# d-subsets).  For each: every residue mod N is a sum of d distinct elements of S.
import itertools, json, time
from fractions import Fraction
seeds = [
  (2, 56,  [3, 4, 8, 10, 12, 16, 17, 28, 31, 38, 45, 48]),
  (2, 6,   [0, 1, 2, 4]),
  (3, 20,  [1, 4, 10, 17, 18, 19]),
  (3, 95,  [7, 10, 26, 40, 64, 65, 75, 83, 85, 90]),
  (3, 63,  [2, 10, 16, 18, 38, 39, 49, 56, 60]),
  (4, 66,  [0, 7, 19, 23, 46, 47, 48, 57]),
  (4, 103, [4, 9, 14, 38, 41, 44, 67, 79, 93]),
  (5, 103, [4, 9, 14, 38, 41, 44, 67, 79, 93]),   # complement of the d=4 set (9 elements, 9-4 = 5)
  (6, 63,  [2, 10, 16, 18, 38, 39, 49, 56, 60]),   # complement of the d=3 set (9 elements, 9-3 = 6)
  (5, 210, [6, 64, 67, 69, 74, 133, 153, 173, 197, 205]),
  (4, 145, [21, 45, 48, 51, 69, 82, 93, 104, 139, 143]),
  (6, 145, [21, 45, 48, 51, 69, 82, 93, 104, 139, 143]),   # complement of the d=4 set (10 elements)
]
out = []
for d, N, S in seeds:
    t0 = time.time()
    ok = len(set(S)) == len(S) and len({sum(c) % N for c in itertools.combinations(S, d)}) == N
    n = len(S)
    out.append(dict(d=d, N=N, n=n, S=S, covers=ok, density=str(Fraction(N, n ** d)), density_float=N / n ** d,
                    runtime_s=round(time.time() - t0, 4)))
    print(d, N, n, ok, Fraction(N, n ** d), round(N / n ** d, 6))
json.dump(dict(command='python3 seeds.py', note='exact check: all d-subsets enumerated with Python ints', seeds=out),
          open('certificates/seeds.json', 'w'), indent=1)
