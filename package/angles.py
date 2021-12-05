import numpy as np
from scipy.optimize import root_scalar
import TaylorMaccoll as tm
import flows as fl

def betaCone(delta: float, gas: fl.Gas)-> float:
    """
    Calc cone shock angle given a semi aperture cone angle

    Input:
    - delta: semi-aperture cone angle in radians
    - gas: gas like dataclass
    - beta_i: initial guess fo shock angle in radians
 
    Returns:
    - beta: cone shock angle in radians
    """
    def obliqueShock(angle):
        x = np.tan(delta)
        for beta in np.arange(1, 500) * np.pi/1000:
            r = 2 / np.tan(beta) * (gas.Ma**2 * np.sin(beta)**2 - 1) / (gas.Ma**2 * (gas.gamma + np.cos(2 * beta)) + 2)
            if r > x:
                break
        return beta
    
    wedgeAngles = np.linspace(delta, 0, 300)
    for wedgeAngle in wedgeAngles:
        beta=obliqueShock(wedgeAngle)
        _,Ma=tm.solveTaylorMaccoll(beta,gas)
        Ma_omega=Ma[1]
        if Ma_omega[-1]<0:
            return beta

def calcMaxDelta(gas: fl.Gas,Mach:np.ndarray = None)-> tuple:
    """
    Calc max semi aperture angle of a cone at zero angle of attack
    and max shock angle given the asintotic Mach

    Parameters:
    - gas: 
    - Mach: optional, numpy ndarray of number of mach at which calculating
            deltaMax and shock angle
    
    Return:
    - delta: ndarray of semi aperture cone angles in radians
    - beta: ndarray of shock angles corresponding to M2=1 in radians
    """

    if isinstance(Mach, type(None)):
        Mach=np.array([gas.Ma])

    # to use also Single Mach values given as float or int
    elif isinstance(Mach,(int,float)):
        Mach=np.array([Mach])
    
    # check is Mach number is subsonic
    if isinstance(Mach,np.ndarray) and Mach.any() <1.0:
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
    beta=np.zeros_like(Mach)

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
            x0=np.pi/3,
            x1=np.pi/4,
            maxiter=100,
            rtol=1e-8)
        
        beta[i]=result.root
        omega,_ = tm.solveTaylorMaccoll(result.root,gas)
        delta[i]=omega[-1]
    
    return delta, beta

if __name__ == "__main__":

    import matplotlib.pyplot as plt
    # dati gas
    Ma =2.1
    dict_air={
    "Ma": Ma,
    "gamma": 1.4,
    "R":287.05,
    "T":273.0,
    "rho":0.129,
    "p": 0.1*101325, 
    "n" : 5}

    air = fl.Gas(dict_air)
    deltac=np.deg2rad(10)
    # beta_0=fl.obliqueShock(deltac,air) ; beta_1=0.8*beta_0
    #beta_0=np.deg2rad(85) ; beta_1=np.deg2rad(60) # with this it converges to the strong solution
    betac=betaCone(deltac,air)

    w,Ma = tm.solveTaylorMaccoll(betac,air)
    Mw=Ma[1]
    Mr=Ma[0]
    print("beta cone = ", np.rad2deg(betac))

    print(f"Ma = {air.Ma}  deltatheta= {np.rad2deg(w[-1]-deltac)}")
    # print(f"beta2D= {np.rad2deg(beta2D)} \n betaCone= {np.rad2deg(beta)}")

    fig, ax1 = plt.subplots()
    ax1.plot(np.rad2deg(w), Mw)
    ax1.set_title(r"grafico $\omega$ - $Ma_{\omega}$")
    ax1.set_ylabel(r"$Ma_{\omega}$")
    ax1.set_xlabel(r"$\omega$")
    ax1.grid()

    fig, ax2 = plt.subplots()
    ax2.plot(np.rad2deg(w), Mr)
    ax2.set_title(r"grafico $\omega$ - $Ma_{r}$")
    ax2.set_ylabel(r"$Ma_{r}$")
    ax2.set_xlabel(r"$\omega$")
    ax2.grid()
    plt.show()

    x=np.linspace(0,np.pi/2,50)
    Mach = 25 - (25 - 1.1)* np.cos(x)

    air=fl.Gas(dict_air)
    delta,_=calcMaxDelta(air)
    delta,beta=calcMaxDelta(air,Mach=Mach)
    delta=np.rad2deg(delta)
    beta=np.rad2deg(beta)

    plt.plot(Mach,delta,"ko",label=r"$\delta$")
    plt.plot(Mach,beta,"bo",label=r"$\beta$")
    plt.legend()
    plt.grid()
    plt.xlabel("Ma")
    plt.show()