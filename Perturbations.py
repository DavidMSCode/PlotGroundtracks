"""Functions for perturbed two body orbits around Earth"""
import numpy as np
from Kepler import rhoe, Re, muEarth
J2 = 1.08263e-3      #Coefficient for second zonal term

def SunSyncInc(a,e=0,unit='radians'):
    #Takes the semi-major axis (and eccentricity if not 0) of an orbit around the Earth and returns the 
    #inclination necessary for a sun-synchronous orbit in degrees.
    cosi = -2*rhoe*a**(7/2)*(1-e**2)**2/(3*J2* Re**2 *muEarth**(1/2))
    i = np.arccos(cosi)
    if unit.lower() == "degrees" or unit.lower() == "degree":
        i = i*180/np.pi
    return i

def SunSyncSMA(i,e=0,unit='radians'):
    #Takes the inclination in radians of an earth orbit and returns the necessary semi-major axis for sun-synchronous orbit.
    if unit.lower() == "degrees" or unit.lower() == "degree":
        i = i*np.pi/180
    a = (-np.cos(i)*((3*J2* Re**2 *muEarth**(1/2)))/(2*rhoe*(1-e**2)**2))**(2/7)
    return a