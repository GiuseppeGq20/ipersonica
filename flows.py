import numpy as np
from scipy.optimize import root_scalar
# try do define an interface for Gas class


class Gas:

    def __init__(self,dictfile: dict) -> None:

        for key in dictfile:
            setattr(self,key,dictfile[key])
        self._setValues()

    def _setValues(self):
        self.cp=self.R* self.gamma /(self.gamma - 1)
        self.a= (self.gamma*self.R*self.T)**0.5 #m/s
        self.H= self.cp*self.T + 0.5*((self.Ma*self.a)**2)
        self.T0=self.H/self.cp

    def update(self):
        self._setValues()

def isentropicFlow(Ma: np.ndarray,gas: Gas):
    """
    calc isentropic relations given initial state of the gas and value of Mach
    number at which evaluate them

    Parameters:
    - Ma: numpy array of Mach number
    - gas: instance of Gas class

    Return:
    - T_T0: temperature ratios array
    - p_p0: pressure ratio array
    - rho_rho0: density ratio array

    Example:
    >>> dict_air={"Ma": 5,"gamma": 1.4,"R":287.05,"T":273.0,"rho":0.129,"p": 0.1*101325,"n" : 5}
    >>> air=Gas(dict_air)
    >>> T_T0,p_p0,rho_rho0=isentropicFlow(np.array([0.8,2]),air)
    >>> np.isclose(T_T0,[0.8865,0.5556],rtol=1e-3).all()==True
        True
    """
    T_T0=(1 + ((gas.gamma -1)/2)*(Ma**2))**-1
    p_p0=(T_T0)**(gas.gamma/(gas.gamma -1))
    rho_rho0=(p_p0)**(1/gas.gamma)
    return T_T0,p_p0,rho_rho0

def obliqueShock(theta: float,gas: Gas)-> float:
    """
    calc 2D schock angle

    Inputs:
    - theta: deviation angle of the flow in radians
    - gas: Gas like object

    Return:
    - beta: schock angle in radians

    Examples:
    shock of a 10 degree step
    >>> dict_air={"Ma": 5,"gamma": 1.4,"R":287.05,"T":273.0,"rho":0.129,"p": 0.1*101325,"n" : 5}
    >>> air1=Gas(dict_air)
    >>> beta=obliqueShock(np.deg2rad(10),air1)
    >>> np.isclose(beta, 0.3382, rtol=1e-3)
        True
    """

    # sanity check
    if gas.Ma <1 : raise RuntimeError("gas is subsonic") 
    # initialization
    elif gas.Ma < 4 : beta0= np.arcsin(1/gas.Ma)
    else : beta0 = theta*((gas.n+1)/(2*gas.n) + np.sqrt( ((gas.n+1)/(2*gas.n))**2 + (gas.Ma*theta)**-2 ))

    beta= root_scalar(
        lambda beta: (1/(1+gas.n))*(1+ gas.n/((gas.Ma**2)*(np.sin(beta)**2))) - np.tan(beta - theta)/np.tan(beta),
        method='secant',
        x0= beta0,
        x1=1.3*beta0,
        maxiter=200,
        xtol=1e-6
    )

    if not beta.converged: raise RuntimeWarning(beta.flag)

    if beta.root > np.deg2rad(90): raise RuntimeError("beta must be less than 90 degrees")
    
    return beta.root

def normalShockRatio(gas: Gas,beta=np.pi/2):
    Ma=gas.Ma*np.sin(beta)
    Ma2 = np.sqrt((Ma**2 + gas.n) / ((gas.n + 2)*(Ma**2) -1 ))
    p2p1= (1 + gas.gamma*(Ma**2)) / (1 + gas.gamma*(Ma2**2))
    rho2rho1=(1/p2p1)* ((Ma/Ma2)**2)
    T2T1=( 1 + (Ma**2)/gas.n) / ( 1 + (Ma2**2)/gas.n)
    
    return Ma2,p2p1,rho2rho1,T2T1

def normalShock(gas: Gas,beta=np.pi/2):
    Ma2,p2p1,rho2rho1,T2T1 = normalShockRatio(gas,beta=beta)
    dict2={
        "Ma": Ma2,
        "gamma": gas.gamma,
        "R":gas.R,
        "T":gas.T*T2T1,
        "rho":gas.rho*rho2rho1,
        "p": gas.p*p2p1, 
        "n" : gas.n
    }
    return Gas(dict2)

if __name__=="__main__":

    import doctest
    doctest.testmod()