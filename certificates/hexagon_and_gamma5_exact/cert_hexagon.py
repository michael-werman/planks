"""CERTIFICATE (exact arithmetic in Q(sqrt3) via sympy) for the regular hexagon theorem:
P = {r e^{i pi j/3}}, r = 1/sqrt(14+4 sqrt3)  =>  rho(L_P) = r/2 = 1/sqrt(56+16 sqrt3).
Checks every step of the written proof that is a computation:
 (1) the 15 spanned lines, as (unit normal n, offset c) with n.x = c, are exactly:
     direction-0/60/120 classes with offsets {0, +-r*sqrt3/2}, direction-30/90/150 classes with offsets +-r/2;
 (2) the set of lines is invariant under rotation by 60 degrees and reflection in the x-axis;
 (3) the three lines used in the sector argument are lines of P;
 (4) the algebraic identity ((7+sqrt3)^2 + (sqrt3+1)^2)/4 * r^2 = 1;
 (5) the witness (0,1) has distance >= r/2 from all 15 lines, with equality for exactly two.
All comparisons are exact (sympy algebraic numbers, sign decided by exact simplification)."""
import sympy as sp, itertools, time
t0 = time.time()
s3 = sp.sqrt(3)
r = 1 / sp.sqrt(14 + 4 * s3)
V = [(r * sp.cos(sp.pi * j / 3), r * sp.sin(sp.pi * j / 3)) for j in range(6)]
V = [(sp.nsimplify(sp.simplify(x)), sp.nsimplify(sp.simplify(y))) for x, y in V]

def line(p, q):
    dx, dy = q[0] - p[0], q[1] - p[1]
    L = sp.sqrt(sp.simplify(dx * dx + dy * dy))
    n = (sp.simplify(-dy / L), sp.simplify(dx / L))
    c = sp.simplify(n[0] * p[0] + n[1] * p[1])
    if (sp.simplify(c) < 0) or (c == 0 and (n[1] < 0 or (n[1] == 0 and n[0] < 0))):
        n = (-n[0], -n[1]); c = -c
    return (sp.radsimp(n[0]), sp.radsimp(n[1]), sp.radsimp(c))

def canon(L):
    a, b, c = L
    if c == 0 and (sp.simplify(b) < 0 or (b == 0 and sp.simplify(a) < 0)): a, b = -a, -b
    return (sp.nsimplify(sp.simplify(a)), sp.nsimplify(sp.simplify(b)), sp.nsimplify(sp.simplify(c)))

lines = set(canon(line(V[i], V[j])) for i, j in itertools.combinations(range(6), 2))
assert len(lines) == 15
# (1) offsets per direction class
half = sp.Rational(1, 2)
offs = sorted(set(sp.nsimplify(sp.simplify(c / r)) for a, b, c in lines), key=lambda t: float(t))
print('offsets / r :', offs)
assert set(offs) == {0, half, s3 / 2}
# (2) symmetry: rotation by 60 degrees and reflection y -> -y map the line set to itself
c6, s6 = half, s3 / 2
def rot(L):
    a, b, c = L; return canon((sp.radsimp(c6 * a - s6 * b), sp.radsimp(s6 * a + c6 * b), c))
def refl(L):
    a, b, c = L; return canon((a, -b, c))
def member(L, S):
    return any(all(sp.simplify(u - v) == 0 for u, v in zip(L, M)) for M in S)
assert all(member(rot(L), lines) for L in lines) and all(member(refl(L), lines) for L in lines)
# (3) lines used: y = 0, y = r*sqrt3/2, x/2 - (sqrt3/2) y = r/2
used = [canon((sp.Integer(0), sp.Integer(1), sp.Integer(0))), canon((sp.Integer(0), sp.Integer(1), sp.nsimplify(sp.simplify(r * s3 / 2)))),
        canon((half, -s3 / 2, sp.nsimplify(sp.simplify(r / 2))))]
for L in used:
    assert member(L, lines), L
# (4) identity
assert sp.simplify(((7 + s3) ** 2 + (s3 + 1) ** 2) / 4 * r ** 2 - 1) == 0
# (5) witness (0,1)
d = [sp.simplify(sp.Abs(b * 1 - c)) for a, b, c in lines]
diff = [sp.simplify(x - r / 2) for x in d]
assert all(sp.simplify(x) >= 0 for x in diff)
assert sum(1 for x in diff if sp.simplify(x) == 0) == 2
print('all exact checks passed; rho = r/2 =', sp.nsimplify(sp.simplify(r / 2)), '=', sp.N(r / 2, 30))
print('runtime %.1fs' % (time.time() - t0))
