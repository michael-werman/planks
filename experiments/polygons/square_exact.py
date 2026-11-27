# Checks (floating point quadrature, mpmath 30 digits) of the closed form for the unit square:
# A(a) = int_0^{pi/2} min(1/sin t, a/cos t)^3 dt; h(x,0) = (2/3)(A(x)+A(1-x));
# claim: min over the boundary is h(1/2,0) = sqrt5/3 + (11/12) log(golden ratio), c = 8/m.
from mpmath import mp, quad, sin, cos, pi, sqrt, log, atan, mpf, diff
mp.dps=30
def A(a):
    if a==0: return mpf(0)
    ts=atan(1/a)
    return quad(lambda t:(a/cos(t))**3,[0,ts])+quad(lambda t:(1/sin(t))**3,[ts,pi/2])
g=(1+sqrt(5))/2
m=sqrt(5)/3+mpf(11)/12*log(g)
print('A(1/2) quad   ', A(mpf(1)/2)); print('A(1/2) closed ', sqrt(5)/4+mpf(11)/16*log(g))
print('m=', m, ' (4/3)A(1/2)=', mpf(4)/3*A(mpf(1)/2), ' c=8/m=', 8/m, ' 96/(4sqrt5+11 log g)=', 96/(4*sqrt(5)+11*log(g)))
for x in ['0','0.1','0.25','0.4','0.5']:
    x=mpf(x); print('h(%s,0)='%x, mpf(2)/3*(A(x)+A(1-x)))
# A'(a) = 3 g(1/a), g(u) = (u sqrt(1+u^2)+asinh u)/(2u^2)
from mpmath import asinh
for a in ['0.2','0.5','0.9']:
    a=mpf(a); u=1/a; print('A\'(%s) num'%a, diff(A,a), ' formula', 3*(u*sqrt(1+u*u)+asinh(u))/(2*u*u))
