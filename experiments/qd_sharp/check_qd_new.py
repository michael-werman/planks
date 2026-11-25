"""Exact (Fractions) checks for the new constants with Lambda'_d(m) = A_d + (1+m)^2 (C_d - 1) + 4(1+m)(1 - bbar).
Criterion: s(1/3 + m^2) - q Lambda'_d(m) is a quadratic in m with positive leading coefficient and negative
discriminant  =>  q_d' >= q/s.  pi < 3.1416, 2.236 <= sqrt5, 1/sqrt6 <= 0.40825, sqrt(pi/2) <= 1.2534."""
from fractions import Fraction as F
import sys
PI_UP = F(31416, 10000)
ok = True
def check(name, c):
    global ok; print(f'{name}: {c}'); ok &= bool(c)
def quad_coeffs(d, Cup, q, s):
    A = F(3*d - 1, d - 1); bb = F(2, d - 1); c = Cup - 1; e = 1 - bb
    # Lambda'(m) = A + c(1+m)^2 + 4e(1+m) = c m^2 + (2c + 4e) m + (A + c + 4e)
    a2, a1, a0 = c, 2*c + 4*e, A + c + 4*e
    Q2, Q1, Q0 = s - q*a2, -q*a1, s/3 - q*a0
    return Q2, Q1, Q0, Q1*Q1 - 4*Q2*Q0
for d, Cup, q, s, nm in [(4, PI_UP/2, F(588, 10000), F(1), "q4' >= 0.0588"), (5, F(2), F(496, 10000), F(1), "q5' >= 0.0496"),
                          (6, 3*PI_UP/4, F(8, 100), F(2), "q6' >= 0.08/2"), (7, F(8, 3), F(8, 100), F(2236, 1000), "q7' >= 0.08/sqrt5")]:
    Q2, Q1, Q0, D = quad_coeffs(d, Cup, q, s)
    print(f'  d={d}: Q2={Q2} ({float(Q2):.6f}) Q1={Q1} ({float(Q1):.6f}) Q0={Q0} ({float(Q0):.6f}) disc={float(D):.6g}')
    check(nm + ' (Q2 > 0, disc < 0)', Q2 > 0 and D < 0)
check('0.0588 >= 0.08/sqrt2  (2*0.0588^2 >= 0.0064)', 2*F(588, 10000)**2 >= F(64, 10000))
check('0.0496 >= 0.08/sqrt3  (3*0.0496^2 >= 0.0064)', 3*F(496, 10000)**2 >= F(64, 10000))
# d >= 8: Lambda'/x <= (23/7) r6 + rp (1+m)^2 + (20/7) r6 (1+m), with C_d - 1 <= C_d <= sqrt(pi/2) x
r6, rp, q = F(40825, 100000), F(12534, 10000), F(8, 100)
Q2 = 1 - q*rp; Q1 = -q*(2*rp + F(20, 7)*r6); Q0 = F(1, 3) - q*(F(23, 7)*r6 + rp + F(20, 7)*r6)
D = Q1*Q1 - 4*Q2*Q0
print(f'  d>=8: Q2={float(Q2):.6f} Q1={float(Q1):.6f} Q0={float(Q0):.6f} disc={float(D):.6f}')
check("d >= 8: q_d' >= 0.08/sqrt(d-2) (Q2 > 0, disc < 0)", Q2 > 0 and D < 0)
check('L4 >= 24*1.0588 = 25.4112, L5 >= 120*1.0496 = 125.952', 24*F(10588, 10000) == F(254112, 10000) and 120*F(10496, 10000) == F(125952, 1000))
sys.exit(0 if ok else 'FAILED')
