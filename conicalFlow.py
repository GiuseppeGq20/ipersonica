# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp, trapezoid
from scipy.optimize import root_scalar
import flows as fl
import copy

def _TaylorMaccoll(w, y, gas):

    v_r, v_w = y  # unpack variables

    A = (gas.gamma - 1)*(1 - v_w**2 - v_r**2)/2
    dydw = [v_w, 
     (v_r*(v_w**2) - A*(2*v_r + (v_w/np.tan(w)))) / (A - v_w**2)]

    return dydw

def solveTaylorMaccoll(beta: float, gas: fl.Gas):
    """
    Solve the Taylor Maccoll equation given a shock wave angle and the type of gas

    Parameters:
    - beta: wave shock angle in radians
    - gas: Gas like structure

    Returns:
    - w: omega angle at wich the solution is calculated
    - Vr,Vw: velocity components,
    """

    #checks
    if beta <0: 
        raise RuntimeError("beta must be positive")

    # calc shock relation needed
    Mn2,p2p1,rho2rho1,T2T1 = fl.normalShockRatio(gas,beta=beta)
    a2=gas.a*np.sqrt(T2T1)
    Vlim=(2*gas.H)**0.5
    Vr_0=gas.Ma*gas.a*np.cos(beta)/Vlim
    Vw_0=-Mn2*a2/Vlim

    # trigger event to end ode integration
    def event(w, y):
        _, v_w = y
        return v_w
    event.terminal = True

    # initial condition
    # note: initial condition are the ones after the shock wave
    y1_0=Vr_0; y2_0=Vw_0
    
    sol = solve_ivp(
        lambda w, y: _TaylorMaccoll(w, y, gas),
        [beta, 0.0],
        [y1_0, y2_0],
        events=event,
        method="RK45",
        rtol=1e-12,
        atol=1e-15
        )

    if not sol.success:
        raise RuntimeError(sol.message) 

    return sol.t ,sol.y*Vlim

def betaCone(delta: float, beta_0: float,beta_1: float, gas: fl.Gas)-> float:
    """
    Calc cone shock angle given a semi aperture cone angle

    Input:
    - delta: semi-aperture cone angle in radians
    - gas: gas like dataclass
    - beta_i: initial guess fo shock angle in radians
 
    Returns:
    - beta: cone shock angle in radians
    """


    def error(beta):
        w,_= solveTaylorMaccoll(beta,gas)
        return w[-1] - delta
    
    # beta_0=fl.betaMax(delta,air)
    # beta_1=0.9*beta_0
    betac= root_scalar(
        error,
        method='secant',
        x0= beta_0, # change it to be beta such that M2 =1, pg.183
        x1= beta_1,
        maxiter=100,
        xtol=1e-10
    )
    
    if betac.root < 0: raise RuntimeError("shock angle can't be negative")
    
    if not betac.converged: raise RuntimeWarning(betac.flag)
    
    return betac.root

def calcMaxDelta(Gas: fl.Gas,Mach:np.ndarray = None)-> tuple:
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
    gas=copy.deepcopy(Gas)

    if isinstance(Mach, type(None)):
        Mach=np.array([gas.Ma])

    # to use also Single Mach values given as float or int
    elif isinstance(Mach,(int,float)):
        Mach=np.array([Mach])
    
    # check is Mach number is subsonic
    if isinstance(Mach,np.ndarray) and Mach.any() <1.0:
        raise RuntimeError("Mach number must be greater than 1")    

    #shock relation used
    def func(beta,gas: fl.Gas):
        Mn2,p2p1,rho2rho1,T2T1=fl.normalShockRatio(gas,beta=beta)
        M2q= Mn2**2 + ((gas.Ma*np.cos(beta))**2)/T2T1
        return M2q - 1

    #initialize array of deltas cone
    delta=np.zeros_like(Mach)
    beta=np.zeros_like(Mach)

    
    for i in range(Mach.size):
        gas.Ma = Mach[i]
        #update properties that depends on Mach
        gas.H= gas.cp*gas.T + 0.5*((gas.Ma*gas.a)**2)
        
        result=root_scalar(
            lambda beta: func(beta,gas),
            method='secant',
            x0=np.pi/3,
            x1=np.pi/4,
            maxiter=100,
            rtol=1e-8)
        
        beta[i]=result.root
        omega,_ = solveTaylorMaccoll(result.root,gas)
        delta[i]=omega[-1]
    
    return delta, beta

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

        return f**(-1.5)
    
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
    Cfx = trapezoid(cp,phi) /(2*np.pi)
    Cfy = trapezoid(cp*np.cos(phi),phi) /((2*np.pi)*np.tan(deltaC))

    # calc aerodynamic coefficients
    Cd= Cfx*np.cos(alpha) + Cfy*np.sin(alpha)
    Cl=-Cfx*np.sin(alpha) + Cfy*np.cos(alpha)
    
    return Cl, Cd

if __name__ == "__main__":

    # dati gas
    Ma = 3
    dict_air={
    "Ma": Ma,
    "gamma": 1.4,
    "R":287.05,
    "T":273.0,
    "rho":0.129,
    "p": 0.1*101325, 
    "n" : 5}

    air = fl.Gas(dict_air)

    beta = np.arcsin(1.2/Ma)

    w,V = solveTaylorMaccoll(beta,air)
    Mw=V[1]
    Mr=V[0]

    deltac=np.deg2rad(15)
    beta_0=fl.obliqueShock(deltac,air) ; beta_1=0.8*beta_0
    # beta_0=np.deg2rad(np.pi/3) ; beta_1=np.deg2rad(np.pi/4) # with this it converges to the strong solution
    betac=betaCone(deltac,beta_0,beta_1,air)

    w,Ma = solveTaylorMaccoll(betac,air)
    Mw=Ma[1]
    Mr=Ma[0]
    print("beta = ", np.rad2deg(betac))

    print(f"Ma = {air.Ma}  theta= {np.rad2deg(w[-1])}")
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

    x=np.linspace(0,np.pi/2,50)
    Mach = 25 - (25 - 1.1)* np.cos(x)

    air=fl.Gas(dict_air)
    delta,_=calcMaxDelta(air)
    delta,beta=calcMaxDelta(air,Mach=Mach)
    delta=np.rad2deg(delta)
    beta=np.rad2deg(beta)

    plt.plot(Mach,delta,"k-",label=r"$\delta$")
    plt.plot(Mach,beta,"b-.",label=r"$\beta$")
    plt.legend()
    plt.grid()
    plt.xlabel("Ma")
    plt.show()
    

# %%
