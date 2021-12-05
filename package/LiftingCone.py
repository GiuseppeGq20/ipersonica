# %%
import numpy as np
import flows as fl
from scipy.integrate import trapezoid
import angles


def deltaLocalCone(deltaC: float, alpha: float, phi: np.ndarray) -> np.ndarray:
    """
    Calc equivalent cone semi aperture for the local cone method

    Parameters:
    - deltaC: cone semiaperture in radians
    - alpha: angle of attack in radians
    - phi: azimutal position on the cone

    Return:
    - deltaEq: equivalent cone semi aperture in radians

    Example:
    we want to evaluate the equivalent semi aperture angle for 10 phi 
    location equidistributed around a cone with deltaC of 10 degrees 
    at 5 degrees of angle of attack
    >>> import numpy as np
    >>> deltaC= np.deg2rad(10); alpha=np.deg2rad(5)
    >>> phi = np.linspace(0, 2*np.pi,10)
    >>> deltaEq= deltaLocalCone(deltaC,alpha,phi)
    """
    deltaEq = np.arcsin(
            np.sin(deltaC)*np.cos(alpha) +
            np.cos(deltaC)*np.sin(alpha)*np.cos(phi)
            )
    return deltaEq

def cpHigh(deltaC: float, alpha:float, Ma: float, phi: np.ndarray) -> np.ndarray:

    """
    Calc cp of the cone using the approximation method of High

    Parameters:
    - deltaC: cone semiaperture in radians
    - alpha: angle of attack in radians
    - Ma: asintotical value of the Mach number
    - phi: azimutal position on the cone

    Return:
    - Cp: value of the Cp at the given phi values

    Example:
    calc High Cp on a cone flying at Mach 3 and 5 degree of angle
    of attack, 
    >>> Ma=3; deltaC= np.deg2rad(10); alpha=np.deg2rad(5)
    >>> phi=np.linspace(0,np.pi,2)
    >>> cpH=cpHigh(deltaC,alpha,Ma, phi=phi)
    """

    # function to calc cp given an arrary of equivalent semi apertures
    def cp(deltaE: np.ndarray)->np.ndarray:
        Cp = 4 * (np.sin(deltaE)**2) * (2.5 + (8/ ((Ma**2 - 1)**0.5) )*np.sin(deltaE)) / (
             1 + (16/ ((Ma**2 - 1)**0.5) )*np.sin(deltaE))
        return Cp
    
    # helper function to tidy the calculation of f Mach_cone 
    def func_MachCone():
        
        delta=deltaLocalCone(deltaC,alpha,phi=np.pi/2)

        if Ma*np.sin(delta) < 1:
            f = Ma*np.cos(delta) * np.sqrt(1- np.sin(delta)/Ma) / np.sqrt(
                (1 + np.exp(-1 - 1.52*Ma*np.sin(delta)) ) *
                (1 + (Ma**2 * np.sin(delta)**2)/4 ) )
        else:
            f = Ma*np.cos(delta) * np.sqrt(1- np.sin(delta)/Ma) / np.sqrt(
                1 + 0.35* ((Ma*np.sin(delta))**1.5))
                    
        f= f**(-1.5)
        return f
    
    # evaluete cp at phi=pi/2
    deltaE_star=deltaLocalCone(deltaC,alpha,phi=np.pi/2)
    Cp_Star=cp(deltaE_star)
    
    #calc cp
    deltaE=deltaLocalCone(deltaC,alpha,phi)
    Cp=cp(deltaE)
    Cp = Cp - (Cp - Cp_Star)*func_MachCone()
    
    return Cp

def calcCLCdCone(deltaC: float, alpha:float , phi:np.ndarray, cp: np.ndarray)->tuple:
    """
    Calc cl and cd of a cone at alpha angle of attack given cp distribution along phi
    hypotesis: cp depends only on azimutal coordinate

    Parameters:
    - deltaC: semi aperture angle of the cone, in radians
    - alpha: angle of attack in radians
    - phi: azimutal coordinates in radians, that span from 0 to 2pi
    - cp: pressure coefficient on the phi coordinates

    Return:
    - cl: lift coefficient
    - cd: drag coefficient

    Example:
    calc cl and cd of a cone using the cp from the High method
    >>> Mach=3; deltaC= np.deg2rad(10); alpha=np.deg2rad(5)
    >>> phi=np.linspace(0,2*np.pi,50)
    >>> cpH=cpHigh(deltaC,alpha,Mach, phi=phi)
    >>> cl,cd= calcCLCdCone(deltaC,alpha,phi,cpH)
    """
    # check if phi span from 0 to 2pi
    if not phi[-1] == 2*np.pi and phi[0]==0.0:
        raise RuntimeError("phi must span from 0 to 2 pi radians")
    
    #check if phi and cp are the same length
    if not phi.size == cp.size:
        raise RuntimeError("phi and cp must have the same length")

    # calc coefficients along coordinate axis
    Cfx = trapezoid(cp,phi) * np.sin(deltaC)/(2*np.pi)
    Cfy = trapezoid(cp*np.cos(phi),phi) * np.cos(deltaC)/(2*np.pi)

    # calc aerodynamic coefficients
    Cd= Cfx*np.cos(alpha) + Cfy*np.sin(alpha)
    Cl=-Cfx*np.sin(alpha) + Cfy*np.cos(alpha)
    
    return Cl, Cd

if __name__ == "__main__":

    import matplotlib.pyplot as plt
    # dati gas
    Ma = 10
    dict_air={
    "Ma": Ma,
    "gamma": 1.4,
    "R":287.05,
    "T":273.0,
    "rho":0.129,
    "p": 0.1*101325, 
    "n" : 5}

    air = fl.Gas(dict_air)  

    # delta equivalent
    deltaC= np.deg2rad(10); alpha=np.deg2rad(5)
    phi = np.linspace(0, 2*np.pi,10)
    deltaEq= deltaLocalCone(deltaC,alpha,phi)

    deltaE=deltaLocalCone(deltaC,alpha,phi=np.pi/2)
    
    #cp High
    Mach=10.0
    deltaC= np.deg2rad(10); alpha=np.deg2rad(5)
    phi=np.linspace(0,2*np.pi,50)
    cpH=cpHigh(deltaC,alpha,Mach, phi=phi)
    plt.plot(np.rad2deg(phi),cpH);  plt.show()

    cl,cd= calcCLCdCone(deltaC,alpha,phi,cpH)
    print(f"cl = {cl}\ncd = {cd}\n")


    

# %%