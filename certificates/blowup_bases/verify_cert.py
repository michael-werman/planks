# Independent exact verifier: reads a certificate JSON (S, N, d) and checks, using only integer
# arithmetic, that S consists of n distinct residues mod N and that EVERY residue mod N is a sum of d
# distinct elements of S (all d-subsets of S are enumerated; nothing from the construction is used).
import sys, json, time
import numpy as np
t0 = time.time()
rec = json.load(open(sys.argv[1]))
S = [int(x) for x in rec['S']]; N = int(rec['N']); d = int(rec['d']); n = len(S)
assert all(0 <= x < N for x in S) and len(set(S)) == n
assert 3 * N < 2 ** 62
A = np.array(S, dtype=np.int64)
hit = np.zeros(N, dtype=np.bool_)
if d == 2:
    for k in range(1, n):
        hit[(A[:k] + A[k]) % N] = True
elif d == 3:
    pairs = np.zeros(0, dtype=np.int64)          # sums of pairs among S[0..k-1]
    for k in range(n):
        if len(pairs): hit[(pairs + A[k]) % N] = True
        pairs = np.concatenate([pairs, (A[:k] + A[k]) % N])
else:
    raise SystemExit('d not supported')
miss = int((~hit).sum())
print(json.dumps(dict(file=sys.argv[1], n=n, N=N, d=d, all_residues_covered=(miss == 0), missing=miss,
                      density=N / n ** d, verify_runtime_s=round(time.time() - t0, 1))))
