import numpy as np
import matplotlib.pyplot as plt
import conicalFlow as cf
import flows as fl

# script to compare wedge and cone shocks

# cone shock layer: M=10, deltac=10 deg
M=10
deltaC=np.deg2rad(10)

dict_air={
"Ma": M,
"gamma": 1.4,
"R":287.05,
"T":273.0,
"rho":0.129,
"p": 0.1*101325, 
"n" : 5}

air=fl.Gas(dict_air)

#calc shock angle
beta_0=fl.obliqueShock(deltaC,air) ; beta_1=0.8*beta_0
beta=cf.betaCone(deltaC,beta_0,beta_1,air)
#solve Taylor Maccol
w,v=cf.solveTaylorMaccoll(beta,air)

T_T2,p_p2,rho_rho2=cf.calcThermoQuantities(v,air) # ratio conical flow

air_after=fl.normalShock(air,beta=beta) # thermo conditions past the shock
T=air_after.T*T_T2
p=air_after.p*p_p2
rho=air_after.rho*rho_rho2

#equivalent wedge
air_after=fl.normalShock(air,beta=beta_0) 
# plotting

#plot v
modv=np.sqrt(v[1]**2 +v[0]**2)
M_wedge=air_after.Ma/np.sin(beta_0-deltaC)
plt.figure() 
plt.plot(np.rad2deg(w),modv,"-.r",label="cone")
plt.hlines(M_wedge*air_after.a,np.rad2deg(w[-1]),1.1*np.rad2deg(w[0]),"k",label="wedge")
plt.ylabel("v[m/s]")
plt.xlabel("w[deg]")
plt.grid()
plt.legend(loc="center right")
plt.show()

#plot Mach number
M=modv/np.sqrt(air.gamma*air.R*T)
plt.figure() 
plt.plot(np.rad2deg(w),M,"-.r",label="cone")
plt.hlines(M_wedge,np.rad2deg(w[-1]),1.1*np.rad2deg(w[0]),"k",label="wedge")
plt.ylabel("M")
plt.xlabel(r"$\omega$[deg]")
plt.grid()
plt.legend(loc="center right")
plt.show()

#plot thermo
thermoq_cone=(T,p/101325,rho)
thermoq_wedge=(air_after.T,air_after.p/101325,air_after.rho)
ylabels=("T[K]","p[atm]",r"$\rho$[Kg/m^3]")

for value_cone,value_wedge,ylabel in zip(thermoq_cone,thermoq_wedge,ylabels):
   plt.figure() 
   plt.plot(np.rad2deg(w),value_cone,"-.r",label="cone")
   plt.hlines(value_wedge,np.rad2deg(w[-1]),1.1*np.rad2deg(w[0]),"k",label="wedge")
   plt.ylabel(ylabel)
   plt.xlabel(r"$\omega$[deg]")
   plt.grid()
   plt.legend(loc="center right")
   plt.show()
