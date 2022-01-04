"""
Script to calculate the relaxation layer downstream a Normal shockWave
using Monti-Napolitano model.

author: Giuseppe Giaquinto
bib: Elementi di Aerodinamica Ipersonica/Rodolfo Monti, Gennaro Zuppardi
"""
# %%
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import flow.flows as fl
import json

# Constants
# units: metric system

#monti napolitano
DELTA=0.1074
BETA=0.75
EPS1=2.891
EPS2=0.2143

#kinetic theory
#eq Tv, nitrogen vibration
C1=7.21434e-4 #[Pa/s]
C2=1.91e6 #[K]
TAU=5.68e-2

#eq alpha, doxigen dissociation
MO2=32 #Kg/Kmole

# ref quantity ,pag.55 cap 3
T_R=59273 #[K]
U_R=3924 #[m/s]
H_R=1.54e7 #[J/Kg]
S_R=260 #[J/Kg K]
P_R=1.45e11 #[Pa]
V_R=1.06e-4 #[m^3/kg]
KR=6.7e8 #[m^6/ kmol^2 s]
T_VC=3364 #[K]

def rhs(x,y,k1):
    """
    define the right hand side of the ODE system to integrate
    """
    u,Tv,alpha,T,rho,p=y

    # evolution vibrational Temperature of N2
    tv=(C1/p)*np.exp((C2/T)**(1/3))

    #dTvdx= (( np.exp(-(T_VC/T) + (T_VC/Tv)) - 1) / (1 - np.exp(-T_VC/T)) ) / (u*tv*TAU)
    dTvdx= ( ( ( (1 - np.exp(-T_VC/Tv)) * np.exp(-T_VC/T + T_VC/Tv) ) / (1 - np.exp(-T_VC/T)) ) -1 )/ (u*tv*TAU)

    # evolution of dissociation of O2
    K= (4/ (MO2*V_R) )* np.exp(- T_R/T) #[kmol/m^3]

    dalphadx= ((2*(1+alpha+DELTA)*(rho**2))* (KR/(MO2**2)) * ( ((1-alpha-BETA)/(4*rho))*K*MO2 - alpha**2 ))/u

    # velocity
    dudx= - U_R*((1 - ((T/T_R)*EPS1)/(1+alpha+DELTA))*dalphadx + EPS2*dTvdx/T_R ) /  ((EPS1 + 1 + alpha + DELTA)*(T/T_R)*((U_R/u) - (k1*V_R*P_R/(p*U_R))) + u/U_R )

    # thermodynamic quantity
    drhodx= - (rho/u) *dudx
    dpdx = - k1*dudx
    dTdx= T*(dudx/u - dalphadx/(1+alpha+DELTA) + dpdx/p)

    dydx=[dudx,dTvdx,dalphadx,dTdx,drhodx,dpdx]
    return dydx

def plot(x,y):
    """
    routine to plot the results
    """
    labels=('u',"Tv","alpha","T","rho","p","alpha_e")
  
    #equilibrium plot,
    x_ep=x[np.argwhere((y[3]-y[1])/ y[1] < 1e-3 )[0]] # point of partial equilibrium
    x_e=x[np.argwhere((y[6]-y[2])/ y[2] < 1e-3 )[0]]
    #temperature and concentrations
    fig,ax1=plt.subplots()
    p1, = ax1.semilogx(x,y[2],color="red", label= r"$\alpha$")
    p2, = ax1.semilogx(x,y[6],color="pink", label= r"$\alpha_e$")
    ax1.grid()

    ax1.set_xlabel("x[m]")

    twin1=ax1.twinx()
    p3, =twin1.semilogx(x,y[1],color="green", label= r"$T_v$")
    p4, =twin1.semilogx(x,y[3],color="blue", label= r"$T$")
    twin1.tick_params(axis='y',width=2,size=4)
    twin1.set_ylabel("[K]")
    ymin,ymax=twin1.get_ylim()
    twin1.vlines(x_ep,ymin,ymax, colors="black", linestyles="dashed") # partial eq
    twin1.vlines(x_e,ymin,ymax, colors="black", linestyles="dashed") # complete eq
    
    ax1.legend(handles=[p1,p2,p3,p4], loc='center left')
    plt.show()
    
    #plot fluiddynamic and thermodynamic properties
    num=(0,3,4,5)
    ylabel=("u [m/s]","T [K]",r"$\rho$ [Kg/m^3]","p [Pa]")
    for i,label in zip(num,ylabel):
        plt.semilogx(x,y[i])
        ymin,ymax=plt.ylim()
        plt.vlines(x_ep,ymin,ymax, colors="black", linestyles="dashed") # partial eq
        plt.vlines(x_e,ymin,ymax, colors="black", linestyles="dashed") # complete eq
        plt.xlabel("x[m]")
        plt.ylabel(label)
        plt.grid()
        plt.show() 



def main():

    # define initial conditions

    filename="45km.json"
    with open(filename) as file:
        gas_dict= json.load(file)

    air=fl.Gas(gas_dict)

    #gas properties downstream of the shock
    air2=fl.normalShock(air)
    u2=air2.Ma*air2.a

    #initial conditions
    y0=[u2, air.T,0, air2.T, air2.rho, air2.p]

    k1=air2.rho*u2

    # integrate ODE system
    result=solve_ivp(
        lambda x,y: rhs(x,y,k1),
       [0,5e4],
       y0,
       method="LSODA",
       rtol=1e-12,
       atol=1e-14
    )

    T=result.y[3]
    rho=result.y[4]
    #alpha equilibrium
    factor=np.exp(-T_R/T) / (rho*V_R)
    alpha_e= -(factor/2) + np.sqrt(factor*(factor/4 + 1 - BETA))

    # plotting
    y=[value for value in result.y]
    y.append(alpha_e)
    plot(result.t,y)


if __name__=="__main__":

    main()

# %%
