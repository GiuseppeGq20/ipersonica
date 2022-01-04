# %%
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import flow.flows as fl
# constants units: metric system

#monti napolitano
DELTA=0.1074
BETA=0.75
EPS1=2.891
EPS2=0.2143

#cinetica chimica
#eq Tv, vibrazione azoto
C1=7.21434e-4 #[Pa/s]
C2=1.91e6 #[K]
TAU=5.68e-2

#eq alpha, dissociazione ossigeno
MO2=32 #Kg/Kmole
# mO2Kg=2.6566962e-26
# k=1.38e-23
R=8.3144
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

def main():
    Ma= 6.4
    gas_dict={
        "Ma": Ma,
        "gamma": 1.4,
        "R":287.05,
        "T":250.0,
        "rho":4e-3,
        "p": 290, 
        "n" : 5
    }
    air=fl.Gas(gas_dict)
    #update gas properties
    air2=fl.normalShock(air)
    u2=air2.Ma*air2.a

    #initial conditions
    y0=[u2, air.T,0, air2.T, air2.rho, air2.p]
    
    print("initial conditions:\n",y0)
    k1=air2.rho*u2

    t_eval=np.linspace(0,5e4,100000)
    result=solve_ivp(
        lambda x,y: rhs(x,y,k1),
       [0,5e4],
       y0,
       method="LSODA" 
    )
    #TODO add evaluation of alpha equilibrium, pg 75 eq (9.8)
    T=result.y[3]
    rho=result.y[4]
    
    #alpha equilibrium
    factor=np.exp(-T_R/T) / (rho*V_R)
    alpha_e= -(factor/2) + np.sqrt(factor*(factor/4 + 1 - BETA))


    fig,ax=plt.subplots()
    labels=("u","Tv","alpha","T","rho","p","alpha_e")
    # ax.set_xscale("log")

    for y in result.y:
        ax.loglog(result.t,y)
    ax.loglog(result.t,alpha_e)
    ax.grid()
    ax.legend(labels)
    plt.show()

    #pressure
    plt.semilogx(result.t, result.y[-1])
    plt.show()
    #densisty
    plt.semilogx(result.t, result.y[4])
    plt.show()
    #temperature
    plt.semilogx(result.t,T)
    plt.semilogx(result.t,result.y[1])
    plt.show()

    #dissociation mass fraction
    plt.semilogx(result.t, result.y[2])
    plt.semilogx(result.t, alpha_e)
    plt.show()

    plt.semilogx(result.t,T)
    plt.show()

    #velocity
    plt.semilogx(result.t,result.y[0])
    plt.show()


if __name__=="__main__":

    main()

# %%
