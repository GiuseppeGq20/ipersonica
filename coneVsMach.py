import matplotlib.pyplot as plt
import numpy as np

import conicalFlow as cf
import flows as fl

"""
script to compare thermo fluid dynamic properties on the supersonic cone
at various free stream Mach numbers
"""

dict_air = {
    "Ma": 2,
    "gamma": 1.4,
    "R": 287.05,
    "T": 273.0,
    "rho": 0.129,
    "p": 0.1*101325,
    "n": 5}

air = fl.Gas(dict_air)

deltac = np.deg2rad([10, 15, 20]).tolist()
Ma = [x for x in range(8, 19,2)]

# calc cp cone at varing delta
for delta in deltac:
    cp = cf.cpCone(Ma, delta, air)
    angle = round(np.rad2deg(delta))
    plt.plot(Ma, cp, "-o", label=f"{angle}°")
    plt.xlabel(r"$M_\infty$")
    plt.ylabel(r"$c_p$")

plt.grid()
plt.legend(loc="center right")
plt.show()

# thermofluiddynamic properties on the cone
T_ratios = []
rho_ratios = []
M_cones = []
for delta in deltac:
    T = []
    rho = []
    Mc = []
    for mach in Ma:
        air.Ma = mach
        air.update()

        beta_0 = fl.obliqueShock(delta, air)
        beta_1 = 0.8*beta_0
        beta = cf.betaCone(delta, beta_0, beta_1, air)

        w, v = cf.solveTaylorMaccoll(beta, air)

        T_T2, p_p2, rho_rho2 = cf.calcThermoQuantities(v, air)

        _, p2_p1, rho2_rho1, T2_T1 = fl.normalShockRatio(air, beta=beta)

        T_T1 = T_T2[-1]*T2_T1
        modv=np.sqrt((v[0]**2 + v[1]**2))
        M= modv[-1]/(air.a * (T_T1)**0.5)
        rho_rho1 = rho_rho2[-1]*rho2_rho1

        T.append(T_T1)
        rho.append(rho_rho1)
        Mc.append(M)

    T_ratios.append(T)
    rho_ratios.append(rho)
    M_cones.append(Mc)

del T, rho, Mc

# density ratio
for delta, rho in zip(deltac, rho_ratios):
    angle = round(np.rad2deg(delta))
    plt.plot(Ma, rho, "-o", label=f"{angle}°")

plt.xlabel(r"$M_\infty$")
plt.ylabel(r"$\frac{\rho}{\rho_\infty}$[kg/m^3]")
plt.grid()
plt.legend()
plt.show()

# Temperature ratios
for delta, T in zip(deltac, T_ratios):
    angle = round(np.rad2deg(delta))
    plt.plot(Ma, T, "-o", label=f"{angle}°")

plt.xlabel(r"$M_\infty$")
plt.ylabel(r"$\frac{T}{T_\infty}$[K]")
plt.grid()
plt.legend()
plt.show()

# cone mach
for delta, Mc in zip(deltac,M_cones):
    angle = round(np.rad2deg(delta))
    plt.plot(Ma, Mc, "-o", label=f"{angle}°")

plt.xlabel(r"$M_\infty$")
plt.ylabel(r"$M_c$")
plt.grid()
plt.legend()
plt.show()

# beta and deltamax at varing mach
x = np.linspace(0, np.pi/2, 50)
Mach = 25 - (25 - 1.1) * np.cos(x)

air = fl.Gas(dict_air)
delta, _ = cf.calcMaxDelta(air)
delta, beta = cf.calcMaxDelta(air, Mach=Mach)
delta = np.rad2deg(delta)
beta = np.rad2deg(beta)

plt.plot(Mach, delta, "k-", label=r"$\delta_{max}$")
plt.plot(Mach, beta, "b-.", label=r"$\beta$")
plt.legend()
plt.grid()
plt.xlabel(r"$M_\infty$")
plt.show()
