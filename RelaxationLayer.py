import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import flow.flows as fl
# constants units: metric system

#monti napolitano
DELTA=0.1074
BETA=0.75
EPS1=2.891
EPS2=0.214

#cinetica chimica
#eq Tv, vibrazione azoto
C1=7.12e-4
C2=1.91e6
TAU=5.68e-2

#eq alpha, dissociazione ossigeno
mO2=32e3 #uma/Kmole
mO2Kg=2.6566962e-26
#k=1.38e-23
R=8.3144
# ref quantity ,pag.55 cap 3
T_R=59273
U_R=3924
H_R=1.54e7
S_R=260
P_R=1.45e11
V_R=1.06e-4
KR=6.7e8 #[m^6/ kmol^2 s]
T_VC=3364

def rhs(x,y,k1):
    
    u,Tv,alpha,T,rho,p=y

    # evolution vibrational Temperature of N2
    tv=(C1/(p*P_R))*np.exp((C2/(T*T_R))**(1/3))
    
    #dTvdx= (( np.exp(-(T_VC/T) + (T_VC/Tv)) - 1) / (1 - np.exp(-T_VC/T)) ) / (u*tv*TAU)
    dTvdx= ( ( ((1 - np.exp(-TAU/Tv) ) * np.exp(-TAU/T)) / ((1 - np.exp(-TAU/T))* np.exp(-TAU/Tv))  ) -1 )/ (u*tv*TAU)

    # evolution of dissociation of O2
    K= (4/ (mO2*V_R) )* np.exp(- 1/T) #[kmol/m^3]

    dalphadx= ((2*(1+alpha+DELTA)*(rho**2)/((mO2*V_R)**2))* KR * ( ((1-alpha-BETA)/(4*rho))*K*mO2*V_R - alpha**2 ))/u

    # velocity
    dudx= - ((1 - (T*EPS1)/(1+alpha+DELTA))*dalphadx + EPS2*dTvdx ) / ( (EPS1 + 1 + alpha + DELTA)*T*((1/u) - (k1/p)) + u )

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
    u2=2040
    air2=fl.normalShock(air)

    y0=[u2/U_R, air.T/T_R,0, air2.T/T_R, air2.rho*V_R, air2.p/P_R]
    k1=air2.rho*u2 * V_R/U_R

    t_eval=np.linspace(0,10000,100000)
    result=solve_ivp(
        lambda x,y: rhs(x,y,k1),
       [0,10e4],
       y0,
       method="RK45" 
    )
    #TODO add evaluation of alpha equilibrium, pg 75 eq (9.8)
    T=result.y[3]
    rho=result.y[4]
    
    #alpha equilibrium
    factor=np.exp(-1/T) / rho
    alpha_e= -factor/2 + np.sqrt(factor*(factor/4 + 1 - BETA))


    fig,ax=plt.subplots()
    labels=("u","Tv","alpha","T","rho","p")
    # ax.set_xscale("log")

    for y in result.y:
        ax.loglog(result.t,y)
    ax.grid()
    ax.legend(labels)
    plt.show()

    plt.semilogx(result.t, result.y[-1]*P_R)
    plt.show()
    plt.semilogx(result.t, result.y[4]/V_R)
    plt.show()
    plt.semilogx(result.t, result.y[2])
    plt.semilogx(result.t, alpha_e)
    plt.show()

if __name__=="__main__":

    main()
