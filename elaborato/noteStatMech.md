# Minimal Notes on statistical mechanics

## A statistical definition of temperature
If a system that is capable of exchanging energy is in *thermal equilibrium* 
with its environment then * it will appear to choose a macroscopic configuration
that maximizes the number of microstates*.
This implicitly assumes that:
- each one of the possible microstates is equally likely to occur
- the microstates of the system are continuallt changing
- given enough time, the system will explore all possible microstates and spend
equal time in each of them (ergodic hypothesis)


This conditions leads to a "statistical" definition of temperature:
$$
    \frac{1}{k_b T}= \frac{d \ln \Omega}{d E}
$$

## ensemlbes
There are three main ensembles that tends to be used in thermal physics:
1. **Microcanonical ensemble** : an ensemble of systems that have each the same 
fixed energy
2. **Canonical ensemble** : an ensebmle of systems each of wich can exchange 
its energy with a large reservoir of heat (Heat bath). This fixes their temperature
3. **Grand canonical ensemble** : an enseble of systems, each of wich can exchange
energy and particles with a large reservoir. This fixes the system's temperature
and chemical potential.

## Canonical ensemble and Boltzmann distribution
Let us consider a system with energy $\epsilon$ in an heat bath. the total Energy 
of system + heat bath is fixed and has value $E$ (thus the system is part of the canonical ensamble with energy $E$). The heat bath has energy $E-\epsilon$.
To each possible energy $\epsilon$ of the system is associated one microstate.
Thus the probability that the system has energy $\epsilon$ is proportional to:
$$
P(\epsilon) \propto \Omega(E-\epsilon)\times 1
$$
we can take the logarithm of the righthand side, and since $\epsilon << E$ we can perfrom a Taylor series expansion around $\epsilon=0$ up to the first order.
$$
\ln(\Omega(E-\epsilon))=\ln(E)-\frac{d\ln(\Omega(E))}{d\epsilon}\epsilon= \ln(\Omega(E)) - \frac{\epsilon}{k_B T}
$$
so that:
$$
\Omega(E-\epsilon)=\Omega(E)e^{-\frac{\epsilon}{k_B T}}
$$
and then we can conclude that:
$$
P(\epsilon) \propto e^{-\frac{\epsilon}{k_B T}}
$$
To get a proper distribution function we have to nomalize it, so the probability of the system to be in a particular microstate $r$ with energy $E_r$ is given by:
$$
P(r)= \frac{e^{-\frac{E_r}{k_B T}}}{Z}
$$
Where $Z = \sum_i e^{-\frac{E_i}{k_B T}}$ is the partition function of the system. It accounts for all the possible energy levels that the system can assume.

>for a similar discussion on the canonical distribution you can also take a look at this [wiki page](https://it.wikipedia.org/wiki/Distribuzione_di_Boltzmann)

## Maxwell-Boltzmann distribution
To get the distribuion function of a system of molecules ,in thermal equlibrium, moving at speed $v_x$ and $v_x + dv_x$ along the $x$ direction we can recall the Boltzmann distribution so that:
$$
P(v_x) \propto e^{\frac{mv_x^2}{2K_BT}}dv_x
$$
since $v_x$ is a continous variable. Hence this relation for the distribution function $g$ is valid:
$$
g(v_x) \propto e^{\frac{mv_x^2}{2K_BT}}
$$
To normalize the distribution so that $\int_{-\infty}^\infty g(v_x)dv_x =1$ we have to evaluate the integral:
$$
\int_{-\infty}^\infty e^{\frac{mv_x^2}{2K_BT}} dv_x= \sqrt{\frac{2\pi k_BT}m}
$$
The expression of the distribution function $g$ is then:
$$
g(v_x)= \sqrt{\frac m {2\pi k_BT}} e^{\frac{mv_x^2}{2K_BT}}
$$
Analogous expression are valid for the distribution function of molecules moving along the other coordinate axis with prescribed velocity.

>inspecting the $g$ distribution expression we recognise that it is a gaussian distribution with mean 0 and variance $\sigma^2= \frac {K_BT} m$

The distribution function $G$ of the fraction of molecules with velocities between ($v_x,v_y,v_z$) and ($v_x + dv_x,v_y + dv_y,v_z + dv_z$) is given by the product of the $g$ distribution functions (we are assuming that v_x,v_y and v_z are stocastically indipendent) so that:
$$
    G(v_x,v_y,v_z)=g(v_x)g(v_y)g(v_z)=\left(\frac m {2\pi k_BT}\right)^\frac 3 2 e^{\frac{m (v_x^2 + v_y^2 + v_z^2)}{2K_BT}}
$$
>Note: this distribution is important when dealing with free and near free molecular flow

> from previous note we recognise that G is the joint distribution of the $g$ distribution functions that are assumed indipendent from each other. this lead to the conclusion that also G is a multivariate gaussian distribution in wich the correlation terms in the variance matrix are equal to zero. 

### speed distribution
Now we ask what is the distribution function $f$ of molecules ,in thermal equilibrium, moving between speed $\mathbf v$ and $\mathbf v + d\mathbf v$, (now we are not fixing each velocity comoponent as in the previous section when discussing about the $G$ distribution function).
The probability of finding a particle in such state is proportional to a Boltmann fatctor (carachterized by the kinetic energy of the particle $\frac 1 2 m v^2$) times the velocity space volume associated with particles at speed between $\mathbf v$ and $\mathbf v + d\mathbf v$. this correspond to a spherical shell of thickness $dv$, $4\pi v^2 dv$.  
Hence:
$$
P(v) \propto v^2  e^{\frac{m (v^2)}{2K_BT}} dv
$$
and $f$:
$$
f(v) \propto v^2  e^{\frac{m (v^2)}{2K_BT}}
$$
to normalize $f$ we have to evaluate this integral:
$$
\int_0^\infty v^2  e^{\frac{m (v^2)}{2K_BT}} dv = \frac 1 4 \sqrt{\frac \pi {(m/2K_BT)^3}}
$$
the final expression for the so called Maxwell distribution is
$$
f(v)= \frac 4 {\sqrt \pi} \left(\frac m {2K_BT}\right)^{3/2} e^{\frac{m (v^2)}{2K_BT}}
$$
>we can see that $G$ and $f$ differs from a multiplicative factor of $\frac 4 {\sqrt \pi}$, and as one might expect for given the probability of finding a particle moving at speed $\mathbf v$ without jointly fixing its component is greater

>The Maxwell distribution solve the boltzmann equation for a system of particles:
>    - in equilibrium
>    - uniformly distributed in space, 
>    - with a fixed or null convective velocity
>    - non interacting particle (i.e. no coulumb forces)
>
>This is also said to be a *hard sphere* approximation of the system of particles

### carachteristic velocities of the Maxwell distribution
- mean velocity $<v> = \sqrt{\frac{8K_BT}{\pi m}}$
- mean squared velocity $<v^2>= \frac {3K_BT}{m}$, this lead to the conclusion that the mean kinetic energy of a particle is $<E_k>=\frac 3 2 K_BT$
- most likely velocity (wich is the velocity at wich the $f$ distribution reach is maximum value) $v_{max} = \sqrt{\frac{2K_BT}m}$

## Mean free path



## NOTE
### **Thermal equilibrium**
- Two or more systems are said to be in thermal equilibrium if their energy content and temperature are no longer changing with times.  
- Systems in thermal equilibrium have the same temperature