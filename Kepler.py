"""Functions for computing various properties of Keplerian Orbits and constants"""
from cmath import nan
import numpy as np
from rotations import *

muEarth = 3.986004418e5     #km^3/s^2 Earth standard gravitational parameter
Re = 6378.1366            #km Earth mean radius
rhoe = 1.99096871e-7 #(rads/s) Mean motion of Earth around Sun

def rv2elm(r,v,mu,tol=1e-10):
    #  Position & Velocity Magnitudes
    R       = np.linalg.norm(r)
    V       = np.linalg.norm(v)

    #  Angular Momentum Vector
    h       = np.cross(r,v);
    H       = np.linalg.norm(h);

    #  Line of Nodes Vector
    nvec    = np.cross([0, 0, 1],h)
    n       = np.linalg.norm(nvec)

    #  Eccentricity Vector
    evec    = 1/mu*((V**2 - mu/R)*r - (np.dot(r,v))*v);
    e       = np.linalg.norm(evec);

    #  Energy
    xi      = (V**2)/2 - mu/R;

    #  Semimajor Axis (a) & Semillatus Rectum (p)
    if abs(1-e) < tol:
        a   = np.inf
        p   = (H**2)/mu
    else:
        a   = -mu/2/xi
        p   = a*(1 - e**2)

    #  Inclination
    i       = np.arccos(h[2]/H)

    # Right Ascension of Ascending Node
    Om      = np.arccos(nvec[0]/n)
    if (nvec[1] < 0):
        Om = 2*np.pi - Om

    # Argument of Perigee
    w = np.arccos((np.dot(nvec,evec))/n/e)
    if (evec[2] < 0):
        w = 2*np.pi - w

    # True Anomaly
    f = np.real(np.arccos(np.dot(evec,r)/R/e))

    if (np.dot(r,v) < 0):
        f = 2*np.pi - f

    #Mean Anomaly & Eccentric Anomaly
    E = 2*np.arctan2(np.sqrt(1-e)*np.tan(f/2),np.sqrt(1+e))
    if E < 0:
        E = 2*np.pi + E
    M = E - e*np.sin(E)
    if M < 0:
        M = 2*np.pi + M

    #  Special Cases
    #  Initialize s
    s1 = nan
    s2 = nan
    s3 = nan
    #  Elliptical Equatorial (ascending node undefined)
    if (i < tol) and (e >= tol):
        s1 = np.arccos(evec[0]/e);
        if (evec[1] < 0):
            s1 = 2*np.pi - s1;   # Longitude of Perigee
        
    # Circular Inclined (perigee undefined)
    elif (i >= tol) and (e < tol):
        s2 = np.arccos(np.dot(nvec,r)/R/n);    # Argument of Latitude
        if (r[2] < 0):
            s2 = 2*np.pi - s2;
        
    # Circular Equatorial (perigee & ascending node undefined)
    elif (i < tol) and (e < tol):
        s3 = np.arccos(r[0]/R);
        if (r[1] < 0):
            s3 = 2*np.pi - s3;    # True Longitude
    return [p, a, e, i, Om, w, f, E, M, [s2,s3,s1]]

def R_ECI2Orbit(Om,i,T):
    R1 = np.array([[np.cos(T),np.sin(T), 0],[-np.sin(T), np.cos(T), 0], [ 0, 0, 1]])
    R2 = np.array([[1, 0, 0],[0, np.cos(i), np.sin(i)],[0, -np.sin(i), np.cos(i)]])
    R3 = np.array([[np.cos(Om), np.sin(Om), 0],[-np.sin(Om), np.cos(Om), 0],[ 0, 0, 1]])

    R = np.matmul(R1,np.matmul(R2,R3))
    return R

def R_Orbit2ECI(Om,i,T):
    R = R_ECI2Orbit(Om,i,T).T
    return R

def elms2rv(p,e,inc,Om,w,f,s=[nan,nan,nan],mu=muEarth,tol=1e-10):
    if e<tol and inc<tol:
        true_long = nan
        #special case circular equatorial
        if not hasattr(s,"__len__"):
            true_long = s
        elif len(s)==3:
            true_long = s[1]
        if np.isnan(true_long):
            raise ValueError("Input elements indicate a circular equatorial orbit, but the true longitude was not provided.")
        else:
            f = true_long
            w = 0
            Om = 0
    elif e<tol and inc>tol:
        #special case circular inclined
        arg_lat = nan
        if not hasattr(s,"__len__"):
            arg_lat = s
        elif len(s)==3:
            arg_lat = s[0]
        if np.isnan(arg_lat):
            raise ValueError("Input elements indicate a circular inclined orbit, but the argument of latitude was not provided.")
        else:
            f = arg_lat
            w = 0
    elif e>tol and inc<tol:
        #special case elliptical equatorial
        true_w = nan
        if not hasattr(s,"__len__"):
            true_w = s
        elif len(s)==3:
            true_w = s[2]
        if np.isnan(true_w):
            raise ValueError("Input elements indicate a elliptical equatorial orbit, but the true longitude of perigee was not provided.")
        else:
            w = true_w
            Om = 0

    #multiply r magnitude with r unit vector in ECI coordinates
    r_pqw = np.array([p*np.cos(f)/(1+e*np.cos(f)),
                    p*np.sin(f)/(1+e*np.cos(f)),
                    0])

    v_pqw = np.array([-(mu/p)**(1/2)*np.sin(f),
                    (mu/p)**(1/2)*(e+np.cos(f)),
                    0])

    R_PQW2IJK = ROT3(-Om)@ROT1(-inc)@ROT3(-w)
    r_vec = R_PQW2IJK@r_pqw.T
    v_vec = R_PQW2IJK@v_pqw.T
  
    return [r_vec, v_vec]


