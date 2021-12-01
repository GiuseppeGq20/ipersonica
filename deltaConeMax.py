"""
script to calc the maximum semiaperture angle of a cone given
an asintotic Mach Number

"""
#%%
import conicalFlow as cf
import flow.flows as fl
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

dict_air = {
    "Ma": 1.0,
    "gamma": 1.4,
    "R": 287.05,
    "T": 273.0,
    "rho": 0.129,
    "p": 0.1*101325,
    "n": 5
}

def calcMaxDelta(Mach:np.ndarray,gas: fl.Gas)-> np.ndarray:
    """
    Calc max semi aperture angle of a cone at zero angle of attack
    given the asintotic Mach
    """

    # to use also Single Mach values given as float or int
    if isinstance(Mach,(int,float)):
        Mach=np.array([Mach])
    
    # check is Mach number is subsonic
    if Mach.any() <1.0:
        raise RuntimeError("Mach number must be greater than 1")
    
    #shock relation used
    def _normalShock(beta: float,gas: object) -> tuple:

        # if beta <0: 
        #     raise RuntimeError("beta must be positive")
            
        n=gas.n
        # normal mach
        Mn=gas.Ma*np.sin(beta)
        # downstream  normal Mach
        Mn2 = ((Mn**2 + n)/((n+2)*(Mn)**2 - 1))**0.5            
        # densisty ratio
        rho1rho2 = (1 + n/(Mn**2))/(n+1)
        # evaluate 2D flow deviation
        theta2D = beta - np.arctan(np.tan(beta)*rho1rho2)
        # evaluate downstream Mach number
        M2 = Mn2/np.sin(beta-theta2D)

        return M2**2

    #initialize array of deltas cone
    delta=np.zeros_like(Mach)
    

    for i in range(Mach.size):

        # gas=Gas(Mach[i])
        # gas=Gas
        gas.Ma = Mach[i]
        #update properties that depends on Mach
        gas.H= gas.cp*gas.T + 0.5*((gas.Ma*gas.a)**2)
        #dictfile["Ma"]=Mach[i]
        # gas=Gas(dictfile)
        
        result=root_scalar(
            lambda beta: _normalShock(beta,gas) - 1,
            method='secant',
            x0=np.pi/2,
            x1=np.pi/4,
            maxiter=100,
            rtol=1e-8)
        
        omega,_ = cf.SolveTaylorMaccoll(result.root,gas)
        delta[i]=omega[-1]
    
    return delta

#------

dict_air={
    "Ma": 1.0,
    "gamma": 1.4,
    "R":287.05,
    "T":273.0,
    "rho":0.129,
    "p": 0.1*101325, 
    "n" : 5
}

# air=fl.Gas(dict_air)
# air=fl.air(1)

#------
x=np.linspace(0,np.pi/2,50)
Mach = 25 - (25 - 1.1)* np.cos(x)

air=fl.Gas(dict_air)
delta=calcMaxDelta(Mach,air)
delta=np.rad2deg(delta)

plt.plot(Mach,delta,"ko")
plt.show()
# %%
