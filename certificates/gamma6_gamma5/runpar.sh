#!/bin/bash
# usage: runpar.sh n THETA_units [K] [d0]  -- runs cert on 4 processes, prints combined summary.
# Referee fix: fails (exit 1, "CERTIFICATE INVALID") unless every rank exits with status 0 and prints
# exactly one DFS summary line; a crashed or killed rank can no longer be hidden in the TOTAL line.
n=$1; T=$2; K=${3:-40}; d0=${4:-14}
s=$(date +%s.%N)
pids=()
for r in 0 1 2 3; do ./cert $n $T $K 4 $r $d0 > cert_n${n}_${T}_r$r.out & pids+=($!); done
fail=0
for r in 0 1 2 3; do
  if ! wait ${pids[$r]}; then echo "rank $r exited with non-zero status" >&2; fail=1; fi
done
e=$(date +%s.%N)
for r in 0 1 2 3; do
  c=$(grep -c '^DFS' cert_n${n}_${T}_r$r.out)
  if [ "$c" != "1" ]; then echo "rank $r printed $c DFS lines (expected 1)" >&2; fail=1; fi
done
cat cert_n${n}_${T}_r*.out
python3 - "$n" "$T" "$s" "$e" <<'PY'
import sys, re, glob
n, T, s, e = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
tot = {}; cpu = 0
for f in sorted(glob.glob(f'cert_n{n}_{T}_r*.out')):
    for line in open(f):
        if line.startswith('DFS') or line.startswith('pre-phase'):
            for k, v in re.findall(r'(\w+)=([\d.]+)', line):
                if k in ('nodes','D','P','J','M','UNRESOLVED'): tot[k] = tot.get(k, 0) + int(v)
                if k == 'cpu': cpu += float(v)
print(f'TOTAL n={n} THETA={T} units = {int(T)/5120:.10f} deg :', tot, f'cpu={cpu:.1f}s wall={e-s:.1f}s')
PY
if [ $fail -ne 0 ]; then echo "CERTIFICATE INVALID: a rank failed" >&2; exit 1; fi
