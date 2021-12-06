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


def normalShockRatio(gas: Gas):
    Ma2 = np.sqrt((gas.Ma**2 + gas.n) / ((gas.n + 2)*(gas.Ma**2) -1 ))
    p2p1= (1 + gas.gamma*(gas.Ma**2)) / (1 + gas.gamma*(Ma2**2))
    rho2rho1=(1/p2p1)* ((gas.Ma/Ma2)**2)
    T2T1=( 1 + (gas.Ma**2)/gas.n) / ( 1 + (Ma2**2)/gas.n)
    
    return Ma2,p2p1,rho2rho1,T2T1

def normalShock(gas: Gas):
    Ma2,p2p1,rho2rho1,T2T1 = normalShockRatio(gas)
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

    # import doctest
    # doctest.testmod()
    Ma=np.linspace(1,20,50)
    dict_air={
    "Ma": 1,
    "gamma": 1.4,
    "R":287.05,
    "T":273.0,
    "rho":0.129,
    "p": 0.1*101325, 
    "n" : 5}
    beta=np.zeros_like(Ma)
    for i in range(Ma.size):
        dict_air["Ma"]=Ma[i]
        air=Gas(dict_air)
        beta[i]=obliqueShock(np.deg2rad(5),air)
    print(beta)
    import matplotlib.pyplot as plt
    plt.plot(beta)
    plt.show()