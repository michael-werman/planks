# Floating point (mpmath); not part of any proof.
# floating point: relative difference between P(dist(z,H)<t)/t and h(z) in the unit disk at t=0.002
from mpmath import mp, quad, mpf, cos, pi, sqrt
mp.dps=30
def F(p): # antiderivative of (1-p^2)^{3/2} on [-1,1], clipped
    p=max(min(p,mpf(1)),mpf(-1)); s=sqrt(1-p*p)
    return (p*s*(5-2*p*p)/8 + mp.asin(p)*3/8)
J=3*pi/8
t=mpf('0.002')
for r in ['0','0.5','0.9','1']:
    r=mpf(r)
    R=quad(lambda phi:(F(r*cos(phi)+t)-F(r*cos(phi)-t))/t,[0,pi/2,pi],maxdegree=10)/(pi*J)
    h=quad(lambda phi:(1-(r*cos(phi))**2)**1.5,[0,pi/2,pi])*2/(pi*J)
    print(r, R, h, (R-h)/h)
