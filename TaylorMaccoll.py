# %%
from inspect import ClassFoundException
import numpy as np
import matplotlib.pyplot as plt
import flow.flows as fl
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar



def _normalShock(beta: float,gas: fl.Gas) -> tuple:
    n=gas.n
    Mn=gas.Ma*np.sin(beta)
    # downstream  normal Mach
    Mn2 = ((Mn**2 + n)/((n+2)*(Mn)**2 - 1))**0.5

    # densisty ratio
    rho1rho2 = (1 + n/(Mn**2))/(n+1)

    # pressure ratio
    p2p1 = ((n+2)*(Mn)**2 - 1)/(n+1)

    # temperature ratio
    T2t1 = (1 + (Mn**2)/n)/(1 + (Mn2**2)/n)

    # evaluate 2D flow deviation
    theta2D = beta - np.arctan(np.tan(beta)*rho1rho2)
    if theta2D < 0:
        raise(RuntimeError("beta > beta_lim"))
    # evaluate downstream Mach number
    M2 = Mn2/np.sin(beta-theta2D)
    # calc downstream sound speed
    a2 = gas.a * np.sqrt(T2t1)

    # calc Vlim
    Vlim = (2*gas.H)**0.5

    # Initial condition for integrating taylor-Maccol eqs., adimensional velocities
    Vw_0 = -M2*a2*np.sin(beta-theta2D)/Vlim  # w component
    Vr_0 = M2*a2*np.cos(beta-theta2D)/Vlim # radial component

    return Mn2, Vw_0, Vr_0 ,Vlim ,a2


def _TaylorMaccoll(w, y, gas):

    v_r, v_w = y  # unpack variables

    A = (gas.gamma - 1)*(1 - v_w**2 - v_r**2)/2
    dydw = [v_w, 
     (v_r*(v_w**2) - A*(2*v_r + (v_w/np.tan(w)))) / (A - v_w**2)]

    return dydw


def SolveTaylorMaccoll(beta: float, gas: fl.Gas):

    # initial condition are the ones after the shock wave

    # calc shock relation
    Mn2, Vw_0, Vr_0 ,Vlim ,a2 = _normalShock(beta,gas)

    # trigger event to end ode integration
    def event(w, y):
        v_r, v_w = y
        return v_w
    event.terminal = True
    event.direction = 0

    #initial condition
    y1_0=Vr_0; y2_0=Vw_0
    
    sol = solve_ivp(
        lambda w, y: _TaylorMaccoll(w, y, gas),
        [beta, 0.0],
        [y1_0, y2_0],
        events=event,
        method="BDF"
        )

    if not sol.success:
        raise(RuntimeError(sol.message))


    return sol.t ,sol.y*Vlim/a2


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
        w,_= SolveTaylorMaccoll(beta,gas)
        return w[-1] - delta
    
    betac= root_scalar(
        error,
        method='secant',
        x0=beta_0,
        x1= beta_1,
        maxiter=200,
        xtol=1e-8
    )

    if not betac.converged: raise RuntimeWarning(betac.flag)
    
    return betac.root

def deltaLocalCone(deltaC: float, alpha: float, phi: np.ndarray) -> np.ndarray:
    """
    Calc euquivalent cone semi aperture for the local cone method

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
    """

    # function to calc cp given an arrary of equivalent semi apertures
    def cp(deltaE: np.ndarray)->np.ndarray:
        Cp = 4 * (np.sin(deltaE)**2) * (2.5 + (8/ (Ma**2 - 1)**0.5 )*np.sin(deltaE)) / (
             1 + (16/ (Ma**2 - 1)**0.5 )*np.sin(deltaE))
        return Cp
    
    # helper function to tidy the calculation of f Mach cone 
    def func_MachCone():

        delta=deltaLocalCone(deltaC,alpha,phi=np.pi/2)

        f = Ma*np.cos(delta) * np.sqrt(1- np.sin(delta)/Ma) * np.power(
                (1 + np.exp(-1 - 1.52*Ma*np.sin(delta)) )*
                (1 + (Ma**2 * np.sin(delta)**2)/2 )
                    ,-0.5)

        f= f**(-1.5)
        return f
    
    # evaluete cp at phi=pi/2
    deltaE_star=deltaLocalCone(deltaC,alpha,phi=np.pi/2)
    Cp_Star=cp(deltaE_star)
    
    #calc cp
    deltaE=deltaLocalCone(deltaC,alpha,phi)
    Cp=cp(deltaE)
    Cp= Cp - (Cp -Cp_Star)*func_MachCone()
    
    return Cp

def calcCLCdCone():

    # return cl,cd
    pass

if __name__ == "__main__":

    # dati gas
    Ma = 3
    air = fl.air(Ma)
    beta = np.deg2rad(30)

    w,Ma = SolveTaylorMaccoll(beta,air)
    Mw=Ma[1]
    Mr=Ma[0]

    deltac=np.deg2rad(20.4)
    beta_0=fl.obliqueShock(deltac,air) ; beta_1=0.98*beta_0
    # beta_0=np.deg2rad(85) ; beta_1=np.deg2rad(86) # with this it converges to the strong solution
    betac=betaCone(deltac,beta_0,beta_1,air)
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
    Ma=3
    phi=np.linspace(0,np.pi,2)
    cpH=cpHigh(deltaC,alpha,Ma, phi=phi)

# %%
