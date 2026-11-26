# Exact (sympy) evaluation of the limit constant c_B for the unit ball B^d (Paper II, ball corollary):
#   h(z) = 2 * avg_{u in S^{d-1}} (1 - <z,u>^2)^k / J(k),  k = (d^2-1)/2,  J(a) = int_{-1}^{1} (1-t^2)^a dt,
#   minimal on the boundary sphere, where h = 2 J(k + (d-3)/2) / (J((d-3)/2) J(k)),
#   c_B = d^2 d! / h_boundary.
import sympy as sp
t=sp.symbols('t')
def J(a): return sp.sqrt(sp.pi)*sp.gamma(a+1)/sp.gamma(a+sp.Rational(3,2))
for d in range(2,7):
    k=sp.Rational(d*d-1,2); a=sp.Rational(d-3,2)
    h0=sp.nsimplify(sp.simplify(2/J(k)))
    hb=sp.simplify(2*J(k+a)/(J(a)*J(k)))
    c=sp.simplify(d**2*sp.factorial(d)/hb)
    print(d, "h(0) =", sp.simplify(h0), "=", sp.N(h0,6), "| h(boundary) =", hb, "=", sp.N(hb,6), "| c_B =", c, "=", sp.N(c,8))
