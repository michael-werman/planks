# usage: python3 make_cert.py <config> t  -> certificates/<config>_t<t>.json (S, N, d, construction parameters)
import sys, json, time, subprocess
import blowup
name, t = sys.argv[1], int(sys.argv[2])
L0, P, d, dirs, c, m0, step = blowup.CONFIGS[name]
t0 = time.time()
C = blowup.build(L0, P, d, dirs, m0 + step * t, c)
rec = dict(config=name, t=t, level_set=L0, P=P, d=d, directions=dirs, offsets=c, m=m0 + step * t,
           a=C['a'], M=C['M'], N=C['N'], T=C['T'], Tbound=C['Tbound'], g=C['g'], n=C['n'],
           density=f"{C['N']}/{C['n']}^{d}", needed_subsets=C['subs'], S=sorted(C['S']),
           build_runtime_s=round(time.time() - t0, 2),
           command=f"python3 make_cert.py {name} {t}")
json.dump(rec, open(f'certificates/{name}_t{t}.json', 'w'))
print(name, t, 'n', C['n'], 'N', C['N'], 'N/n^d', C['N'] / C['n'] ** d)
