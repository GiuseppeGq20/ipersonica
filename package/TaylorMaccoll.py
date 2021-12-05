import numpy as np
from scipy.integrate import solve_ivp 
import flows as fl

def _normalShock(beta: float,gas: fl.Gas) -> tuple:

    if beta <0: 
        raise RuntimeError("beta must be positive")
     
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

    # if theta2D < 0:
    #     raise(RuntimeError("beta > beta_lim"))
    # evaluate downstream Mach number
    Mt2= gas.Ma*np.cos(beta)/np.sqrt(T2t1)
    M2= (Mn2**2 + Mt2**2)**0.5
    # M2 = Mn2/np.sin(beta-theta2D)
    # calc downstream sound speed
    a2 = float(gas.a) * np.sqrt(T2t1)

    # calc Vlim
    Vlim = (2*gas.H)**0.5

    # Initial condition for integrating taylor-Maccol eqs., adimensional velocities
    Vw_0 = -Mn2*a2/Vlim  # w component
    Vr_0 = Mt2*a2/Vlim # radial component

    return Mn2, Vw_0, Vr_0 ,Vlim ,a2

def taylorMaccoll(w, y, gas):

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

    # calc shock relation
    _, Vw_0, Vr_0 ,Vlim ,a2 = _normalShock(beta,gas)

    # trigger event to end ode integration
    def event(w, y):
        _, v_w = y
        return v_w
    event.terminal = True

    # initial condition
    # note: initial condition are the ones after the shock wave
    y1_0=Vr_0; y2_0=Vw_0
    
    sol = solve_ivp(
        lambda w, y: taylorMaccoll(w, y, gas),
        [beta, 0.0],
        [y1_0, y2_0],
        events=event,
        method="BDF"
        ) #BDF backward difference to have a better resolution

    if not sol.success:
        raise RuntimeError(sol.message) 

    return sol.t ,sol.y*Vlim/a2

if __name__ == "__main__":
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

    beta = np.arcsin(1.2/Ma)

    w,Ma = solveTaylorMaccoll(beta,air)
    Mw=Ma[1]
    Mr=Ma[0]