import numpy as np
from scipy.integrate import solve_ivp
import flow.flows as fl
# constants units: metric system

#monti napolitano
DELTA=0.1074
BETA=0.75
EPS1=2.891
EPS2=0.213

#cinetica chimica
#eq Tv, vibrazione azoto
C1=7.21e-4
C2=1.91e6
TAU=5.68e-2

#eq alpha, dissociazione ossigeno
mO2=2,6566962e-26
k=1.38e-23
R=8.3144
# ref quantity ,pag.55 cap 3
T_R=59273
U_R=3924
H_R=1.54e7
S_R=260
P_R=1.45e11
V_R=1.06e-4

T_VC=3364


def rhs(x,y,k1):
    
    u,Tv,alpha,T,rho,p=y

    # evolution vibrational Temperature of N2
    tv=(C1/p)*np.exp((C2/T)**(1/3))
    
    dTvdx= (( np.exp(-(T_VC/T) + (T_VC/Tv)) - 1) / (1 - np.exp(-T_VC/T)) ) / (u*tv*TAU)

    # evolution of dissociation of O2
    K= (4/ (mO2*V_R))* np.exp(- T_R/T)
    
    KR=A* np.exp(-E_A/(k*T))

    dalphadx= (2*(1+alpha+DELTA)*(rho**2)/(mO2**2))*KR*( ((1-alpha-BETA)/(4*rho))*K*mO2 -alpha**2 )/u

    # velocity
    dudx= - ((1 - (T*EPS1)/(1+alpha+DELTA))*dalphadx + EPS2*dTvdx ) / ( (EPS1 + 1 + alpha + DELTA)*T*((1/u) - (k1/p)) + u )

    # thermodynamic quantity
    drhodx= - (rho/u) *dudx
    dpdx = - k*dudx
    dTdx= T*(dudx/u - dalphadx/(1+alpha+DELTA) + dpdx/p)
    
    dydx=[dudx,dTvdx,dalphadx,dTdx,drhodx,dpdx]
    return dydx

#TODO implement custom solver like in the book
def solve():
    while 

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
    k1=air2.rho*u2

    y0=[u2,air.T,0,air2.T,air2.T,air2.p]

    result=solve_ivp(
        lambda x,y: rhs(x,y,k1),

    )


if __name__=="__main__":
    main()