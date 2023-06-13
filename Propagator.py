from scipy.integrate import solve_ivp
import numpy as np

from Kepler import Re

def J2Accel(t,y):
    J2 = 1.0826e-3
    Req = 6378
    mu = 3.986004418e5     #km^3/s^2 Earth standard gravitational parameter
    [x,y,z,u,v,w] = y

    r = (x**2+y**2+z**2)**(1/2)
    J2term = 3/2*J2*(Req/r)**2
    z2term = 5*z**2/r**2
    muterm = mu/r**3
    
    xdot = u 
    ydot = v
    zdot = w
    udot = -muterm*x*(1+J2term*(1-z2term))
    vdot = -muterm*y*(1+J2term*(1-z2term))
    wdot = -muterm*z*(1+J2term*(3-z2term))

    return [xdot,ydot,zdot,udot,vdot,wdot]

def twoBody(t,y):
    mu = 3.986004418e5     #km^3/s^2 Earth standard gravitational parameter
    [x,y,z,u,v,w] = y
    r = (x**2+y**2+z**2)**(1/2)
    muterm = mu/r**3
    xdot = u 
    ydot = v
    zdot = w
    udot = -muterm*x
    vdot = -muterm*y
    wdot = -muterm*z
    return [xdot,ydot,zdot,udot,vdot,wdot]

def J2Propagate(r0,v0,t0,tf,dt=30):
    y0 = [r0[0],r0[1],r0[2],v0[0],v0[1],v0[2]]
    ts = np.arange(t0,tf,dt)
    sol = solve_ivp(J2Accel,[t0,tf],y0,t_eval=ts,atol=1e-12,rtol=1e-12)
    return sol

def twoBodyPropagate(r0,v0,t0,tf,dt=30):
    y0 = [r0[0],r0[1],r0[2],v0[0],v0[1],v0[2]]
    ts = np.arange(t0,tf,dt)
    sol = solve_ivp(twoBody,[t0,tf],y0,t_eval=ts,atol=1e-12,rtol=1e-12)
    return sol