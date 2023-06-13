"""Functions for plotting the groundtrack of a satellite around the Earth"""
import spiceypy as sp
import numpy as np
Re = 6378.1366            #km Earth mean radius

def SplitGroundtrack(ts,lats,lons,tolerance=180):
    #expand coords to lists of lon and lat
    T=ts
    #find every large jump in value (i.e. the satellite crossed the anti-meridian)
    indx=[]
    for i,x in enumerate(lons):
        if i>0:
            if abs(lons[i]-lons[i-1])>=tolerance:
                #large jump occured in longitude add index
                indx.append(i)
    if len(indx)>0:
        #Satellite crossed the anti meridian at least once so we need seperate list for the groundtracks before and after each crossing
        indx.insert(0,0)        #prepend the first index groundtrack
        indx.append(len(lons))     #append the last index of the groundtack
        split_Lon=[]
        split_Lat=[]
        split_T=[]
        for j,id1 in enumerate(indx):
            if id1>0:
                id0=indx[j-1]
                split_Lon.append(lons[id0:id1])  #split lat lon and time lists into each rev of the globe
                split_Lat.append(lats[id0:id1])
                split_T.append(T[id0:id1])
    else:
        split_Lon=[lons]
        split_Lat=[lats]
        split_T=[T]
    return [split_T,split_Lon,split_Lat]

def wrapToPi(lam):
    #wraps radian value such that the values range from -pi to pi
    if hasattr(lam,'__len__'):
        mask = (lam > np.pi) + (lam < -np.pi)
        lam[mask] = (lam[mask] + np.pi) % (2 * np.pi) - np.pi
    else:
        if (lam > np.pi) or (lam < -np.pi):
            lam = (lam + np.pi) % (2 * np.pi) - np.pi
    return lam

def ECI2ECEF(t,x,y,z):
    """Convert from inertial earth centered coordinates to earth centered earth fixed coordinates"""
    A = sp.pxform("J2000","IAU_EARTH",t)
    pos_ecef = A@np.array([x,y,z])
    return pos_ecef

def ECIs2LonLats(T,X,Y,Z):
    #converts vectors of eci positions and times to latitudes and longitudes for the Earth
    f = 1/298.25700617731906        #Earth flattening coefficient
    lats = []
    lons = []
    lon_rots = []
    for t,x,y,z in zip(T,X,Y,Z):
        pos_ecef = ECI2ECEF(t,x,y,z)
        lon,lat,_ = sp.recgeo(pos_ecef,Re,f)
        lats.append(lat*180/np.pi)
        lons.append(lon*180/np.pi)
    return [lons,lats]