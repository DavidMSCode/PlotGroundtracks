import numpy as np
def ROT1(a,unit="radians"):
    #rotation around 1st principal axis
    if unit.lower()=="degree" or unit.lower()=="degrees":
        #convert input to radians if provided in degrees
        a = a*np.pi/180
    R1 = np.array([[1, 0, 0],[0, np.cos(a), np.sin(a)],[0, -np.sin(a), np.cos(a)]])
    return R1

def ROT2(a,unit="radians"):
    #rotation around 2nd principal axis
    if unit.lower()=="degree" or unit.lower()=="degrees":
        #convert input to radians if provided in degrees
        a = a*np.pi/180
    R2 = np.array([[np.cos(a), 0, -np.sin(a)],[0, 1, 0],[np.sin(a), 0, np.cos(a)]])
    return R2

def ROT3(a,unit="radians"):
    #rotation around 3rd principal axis
    if unit.lower()=="degree" or unit.lower()=="degrees":
        #convert input to radians if provided in degrees
        a = a*np.pi/180
    R3 = np.array([[np.cos(a), np.sin(a), 0],[-np.sin(a), np.cos(a), 0],[0, 0, 1]])
    return R3

    