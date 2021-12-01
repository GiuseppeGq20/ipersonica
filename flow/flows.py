import numpy as np
from scipy.optimize import root_scalar
from dataclasses import InitVar, dataclass, field


@dataclass
class GasDataclass:
    Ma:float
    gamma: float = field(init=False)
    R: float = field(init=False) #Nm/kgK
    T: float = field(init=False) # K
    rho: float = field(init=False) #kg/m3
    p: float = field(init=False) #Pa 
    n: int = field(init=False) #dof
    cp: float = field(init=False)
    a: float = field(init=False)#m/s
    H: float = field(init=False) #J

    def __post_init__(self):
        self.cp=self.R* self.gamma /(self.gamma - 1)
        self.a= (self.gamma*self.R*self.T)**0.5 #m/s
        self.H= self.cp*self.T + 0.5*((self.Ma*self.a)**2) #J

class Gas:

    def __init__(self,dictfile: dict) -> None:

        for key in dictfile:
            setattr(self,key,dictfile[key])
        
        self.cp=self.R* self.gamma /(self.gamma - 1)
        self.a= (self.gamma*self.R*self.T)**0.5 #m/s
        self.H= self.cp*self.T + 0.5*((self.Ma*self.a)**2)



def isentropicFlow():
    pass

def obliqueShock(theta: float,gas: Gas)-> float:
    """
    calc 2D schock angle

    Inputs:
    - theta: deviation angle of the flow in radians
    - gas: Gas like dataclass

    Return:
    - beta: schock angle in radians

    Examples:
    shock of a 10 degree step
    >>> air1=air(3)
    >>> beta=obliqueShock(0.17,air1); round(beta,4)
        0.4779

    """

    # sanity check
    if gas.Ma <1 : raise RuntimeError("gas is subsonic") 
    # initialization
    elif gas.Ma < 4 : beta0= np.arcsin(1/gas.Ma)
    else : beta0 = theta*((gas.n+1)/(2*gas.n) + np.sqrt( ((gas.n+1)/(2*gas.n))**2 + (gas.Ma*theta)**-2 ))

    beta= root_scalar(
        lambda beta: (1/(1+gas.n))*(1+ gas.n/((gas.Ma**2)*np.sin(beta)*np.sin(beta))) - np.tan(beta - theta)/np.tan(beta),
        method='secant',
        x0= beta0,
        x1=1.3*beta0,
        maxiter=200,
        xtol=1e-4
    )

    if not beta.converged: raise RuntimeWarning(beta.flag)
    
    return beta.root




if __name__=="__main__":

    import doctest
    doctest.testmod()
    