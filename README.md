# Planks through points — companion repository

Companion to the paper *Planks through points: covering a ball by slabs around spanned hyperplanes*.

## Contents

- `paper/` — LaTeX sources (pdfLaTeX, main document `main.tex`) and figures.
- `certificates/gamma6_gamma5/` — the only computer-assisted results of the paper: the certified bounds on
  γ₅ and γ₆ (Table 1). Exact integer branch-and-bound for the lower bounds (`cert.c`), interval arithmetic for
  the upper bounds (`upper.py`). `RECORD.txt` gives the commands, parameters, runtimes, machine and outputs;
  `RERUN.txt` a rerun of the γ₆ ≥ 12.55° run after the final code fixes (identical node counts); `sha256.txt` the file hashes. The integer angle table
  (`table_921600.bin`, about 15 MB) is not stored: regenerate it with `python3 mktable.py 921600` and verify it
  with `checktable.py` as described in `RECORD.txt`.
- Independent exact re-checks (Python integers / fractions / sympy; not needed for the proofs, which are by hand):
  - `certificates/hexagon_and_gamma5_exact/` — the hexagon value of Theorem D (ρ² exact in Q(√3)), and γ₅.
  - `certificates/blowup_bases/` — the seeds and two blown-up instances of the restricted additive bases.
  - `certificates/sigma2_seven_eighteenths/` and `proofs/sigma2_seven_eighteenths.md` — the seed behind
    σ₂ ≥ 7/18 (Appendix B.3).
  - `experiments/planar_hand/`, `experiments/theoremA_d3_hand/`, `experiments/qd_sharp/check_qd_new.py` — the
    algebra of Appendix A (constants 6864/3161, 20/3, q_d ≥ 0.08/√(d−2), L₄, L₅).
- `experiments/flats_middle_dims.py` — floating-point exploration for Appendix D.3 (not part of any proof).

## Rerunning the exact checks

```
cd certificates/sigma2_seven_eighteenths && python3 sigma2_check.py      # about 80 s
cd experiments/qd_sharp && python3 check_qd_new.py
cd experiments/planar_hand && python3 verify_identities.py
cd experiments/theoremA_d3_hand && python3 verify_identities.py
cd certificates/hexagon_and_gamma5_exact && python3 cert_hexagon.py       # about 25 s
```

Requirements: Python 3 with sympy, numpy, mpmath (and scipy for the floating-point explorations); gcc for
`cert.c`.
