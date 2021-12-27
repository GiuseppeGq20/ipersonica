---
title: "Minimal Notes on statistical mechanics"  
author: "Giuseppe Giaquinto"
---

# Statistical mechanics

## A statistical definition of temperature

If a system that is capable of exchanging energy is in *thermal equilibrium* 
with its environment then *it will appear to choose a macroscopic configuration
that maximizes the number of microstates*.
This implicitly assumes that:

- each one of the possible microstates is equally likely to occur
- the microstates of the system are continuallt changing
- given enough time, the system will explore all possible microstates and spend
  equal time in each of them (ergodic hypothesis).

This conditions leads to a "statistical" definition of temperature:  
$$
    \frac{1}{k_b T}= \frac{d \ln \Omega}{d E}
$$

## Ensembles

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
of system + heat bath is fixed and has value $E$ (thus the system is part of the 
canonical ensamble with energy $E$). The heat bath has energy $E-\epsilon$.
To each possible energy $\epsilon$ of the system is associated one microstate.
Thus the probability that the system has energy $\epsilon$ is proportional to:
$$
P(\epsilon) \propto \Omega(E-\epsilon)\times 1
$$
we can take the logarithm of the righthand side, and since $\epsilon << E$ we can 
perfrom a Taylor series expansion around $\epsilon=0$ up to the first order.
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
To get a proper distribution function we have to nomalize it, so the probability
of the system to be in a particular microstate $r$ with energy $E_r$ is given by:
$$
P(r)= \frac{e^{-\frac{E_r}{k_B T}}}{Z}
$$
Where $Z = \sum_i e^{-\frac{E_i}{k_B T}}$ is the partition function of the system. 
It accounts for all the possible energy levels that the system can assume.

> for a similar discussion on the canonical distribution you can also take a look 
> at this [wiki page](https://it.wikipedia.org/wiki/Distribuzione_di_Boltzmann)

## Maxwell-Boltzmann distribution

To get the distribuion function of a system of molecules ,in thermal equlibrium, 
moving at speed $v_x$ and $v_x + dv_x$ along the $x$ direction we can recall the 
Boltzmann distribution so that:
$$
P(v_x) \propto e^{\frac{mv_x^2}{2K_BT}}dv_x
$$
since $v_x$ is a continous variable. Hence this relation for the distribution 
function $g$ is valid:
$$
g(v_x) \propto e^{\frac{mv_x^2}{2K_BT}}
$$
To normalize the distribution so that $\int_{-\infty}^\infty g(v_x)dv_x =1$ we 
have to evaluate the integral:
$$
\int_{-\infty}^\infty e^{\frac{mv_x^2}{2K_BT}} dv_x= \sqrt{\frac{2\pi k_BT}m}
$$
The expression of the distribution function $g$ is then:
$$
g(v_x)= \sqrt{\frac m {2\pi k_BT}} e^{\frac{mv_x^2}{2K_BT}}
$$
Analogous expression are valid for the distribution function of molecules moving 
along the other coordinate axis with prescribed velocity.

> nspecting the $g$ distribution expression we recognise that it is a gaussian 
> distribution with mean 0 and variance $\sigma^2= \frac {K_BT} m$

The distribution function $G$ of the fraction of molecules with velocities between
($v_x,v_y,v_z$) and ($v_x + dv_x,v_y + dv_y,v_z + dv_z$) is given by the product
of the $g$ distribution functions (we are assuming that v_x,v_y and v_z are stocastically 
indipendent) so that:  
$$
    G(v_x,v_y,v_z)=g(v_x)g(v_y)g(v_z)=\left(\frac m {2\pi k_BT}\right)^\frac 3 2 e^{\frac{m (v_x^2 + v_y^2 + v_z^2)}{2K_BT}}
$$

> Note: this distribution is important when dealing with free and near free molecular flow

> from previous note we recognise that G is the joint distribution of the $g$ 
> distribution functions that are assumed indipendent from each other. this lead 
> to the conclusion that also G is a multivariate gaussian distribution in wich 
> the correlation terms in the variance matrix are equal to zero. 

### speed distribution

Now we ask what is the distribution function $f$ of molecules ,in thermal equilibrium, 
moving between speed $\mathbf v$ and $\mathbf v + d\mathbf v$, (now we are not 
fixing each velocity component as in the previous section when discussing about 
the $G$ distribution function).
The probability of finding a particle in such state is proportional to a Boltzmann
factor (characterized by the kinetic energy of the particle $\frac 1 2 m v^2$) 
times the velocity space volume associated with particles at speed between
$\mathbf v$ and $\mathbf v + d\mathbf v$. this correspond to a spherical shell 
of thickness $dv$, $4\pi v^2 dv$.  
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

> we can see that $G$ and $f$ differs from a multiplicative factor of $\frac 4 {\sqrt \pi}$, 
> and as one might expect the probability of finding a particle moving at speed 
> $\mathbf v$ without jointly fixing its component is greater

> The Maxwell distribution solve the Boltzmann equation for a system of particles:
> 
> - in equilibrium
> - uniformly distributed in space, 
> - with a fixed or null convective velocity
> - non interacting particle (i.e. no Coulomb forces)
> 
> This is also said to be a *hard sphere* approximation of the system of particles

### characteristic velocities of the Maxwell distribution

we can define the followings:

- mean velocity  $\langle v \rangle = \sqrt{\frac{8K_BT}{\pi m}}$
- mean squared velocity $\langle v^2 \rangle= \frac {3K_BT}{m}$, this lead to the conclusion 
  that the mean kinetic energy of a particle is $\langle E_k\rangle=\frac 3 2 K_BT$
- most likely velocity (wich is the velocity at wich the $f$ distribution reach 
  is maximum value)  $v_{max} = \sqrt{\frac{2K_BT}m}$

## Mean Collision time

the events consistiong of particle colliding into each other distribute as an 
exponential random variate of parameter $n\sigma v$ where $\sigma$ is the 
collisional cross section.Thus the mean collision time is:
$$
\tau= \frac 1 {n\sigma v}
$$

### proof with scattering theory

Let us define $P(t)$ as follows:
$$ P(t)= \text{the probability of not colliding up to time}\quad t $$
In a time $dt$ a molecule will sweep out a volume $\sigma dt$. With $n$ molecules
per unit volume, the probability of a collision in time $dt$ is therefore $n\sigma v dt$.
Where $v $ is the relative velocity between molecules.
Now we can derive the probability of not colliding up to time $t$. We begin by
deriving $P(t + dt)$ wich is the probability of not colliding up to time $t+dt$
$$
P(t+dt)=P(t)P(dt)=P(t)(1-n\sigma v dt)\\
$$ 
but also
$$
P(t+dt)=P(t) + \frac{dP}{dt}dt
$$
this leads us to the following differential equation, with initial condition $P(0)=1$,
$$ \frac{dP}{dt}=-Pn\sigma v$$
and thus:
$$ P(t)= e^{-n\sigma v t}$$
Now the probability of surviving without collision up to time $t$ but then colliding
in the next dt is:
$$ P(t)= e^{-n\sigma v t}n\sigma v dt$$
This define the collision time, wich distribute as an exponential random variate.
We can calculate:

- mean collisional time: $\tau=\frac 1 {n\sigma v}$
- root mean square collisional time: 
  $$\tau_{rms}=\sqrt{\langle t^2\rangle}=(Var(t)+\tau^2)^{1/2}=\frac {\sqrt 2}{n\sigma v}$$

### collisional cross section

We consider particles subjected to a **hard sphere potential** 
$$
V(R)=
\begin{cases}
0 \qquad R >b\\
\infty \qquad R\leq b
\end{cases}
$$
wich implies that two molecule collide if their distance is less than  the impact 
parameter $b$.  
Thus molecules can be tought of moving from one collision to another inside an 
immaginary tube of cross-section $\sigma = \pi b^2$.

> the hard sphere potential approximation is valid for not too low, nor to high 
> temperature

## Mean free path

Consider a first class of molecules wich moves at velocity $\mathbf v$ and
consider only collision with a second class of molecules which move at velocity 
$\mathbf u$ . In a frame moving at velocity $\mathbf u$ , this second class of 
molecules are stationary and offer a total cross-section of 
$n\sigma f(\mathbf u)d\mathbf u$, where $f(\mathbf u)=g(u_x)g(u_y)g(u_z)$ is 
a Maxwell-Boltzmann distribution for the vector $\mathbf u=(u_x,u_y,u_z)$. 
In unit time, the total volume swept out by these targets relative to the 
first class of molecules (which in this frame move at velocity 
$\mathbf v - \mathbf u$) is $| \mathbf v - \mathbf u|n\sigma f(\mathbf u)d\mathbf u$.  
The number of encounters per second is obtained by multipliyng this volume by the 
probability of finding one of the first class of molecules in unit volume,
giving $| \mathbf v - \mathbf u|n\sigma f(\mathbf u)d\mathbf u f(\mathbf v)d\mathbf v$.
The collision rate $R$ is therefore obtained by integrating over all $\mathbf u$
and $\mathbf v$ giving:  

$$ 
R= n \sigma \int\int | \mathbf v - \mathbf u|n\sigma f(\mathbf u)d\mathbf u f(\mathbf v)d\mathbf v
$$  

which writing 
$x=\frac {\mathbf v - \mathbf u}{\sqrt 2}$ and $y=\frac {\mathbf v + \mathbf u}{\sqrt 2}$ 
can be trasformed into  

$$ 
R= n \sigma \sqrt 2 \int |x|n\sigma f(\mathbf x)dx \int f(\mathbf y)d\mathbf y 
$$  

where the first integral yields $\langle v\rangle$ and  the second integral is unity. Hence
$R=n\sigma\sqrt 2 \langle v\rangle$ and the mean free path is:
$$
\lambda = \frac {\langle v\rangle} R = \frac 1 {n \sigma \sqrt 2} 
$$

# Entropy

## thermodynamic definition of entropy

entropy is defined as:
$$
dS = \frac {dQ_{rev}}{T}
$$
and it is a function of state.

> note: since **entropy** is a function of state, its evaluation doesn't depend
> on the particular trasformation to which the system is subjected to. Hence even 
> for an irreversible process we can evaluate the change in entropy of the process
> from state **A** to state **B** by evaluating a the change in entropy of a
> reversible trasformation of the same system that has the same intial and final 
> states **A** and **B**.

### another statement of the II law using the entropy

For a thermally isolated system  the change in entropy is always greater or equal to zero.
$$ dS \geq 0 $$

## Statistical definition of entropy

### Boltzmann definition

if we assume that the statistical definition of temeperature is valid then:
$$
    \frac{1}{k_b T}= \frac{d \ln \Omega}{d E}
$$
and also, from the first law:
$$
\frac 1 T = \left(\frac {\partial S}{\partial U}\right)_V
$$
by comparing these two equation we get:
$$
    S=k_B\ln(\Omega)
$$

> this is a very profound equation, it express the Entropy of a system in one
> possible Macrostate characterized by $\Omega$ possible microstate.

However this isn't the most general case, since we are not accounting the fact
that some system can have more than one possible macrostate (microstate subgroups
that we can address by a macroscopic quantity  or condition of the system) and 
for each macrostate there are $\Omega_i$ microstate, which often are difficult
to measure.  
So in most cases the Boltzmann formula for the entropy would read:
$$
    S_{tot}=k_B\ln(N)
$$  
where $N$ is the number of total microstates. This formula isn't of practical use
due to the difficulties of counting all the microstates.

### Gibbs formula

To overcome this difficulties of using the Boltzmann entropy relation we can use
the Gibbs formula that gives the entropy that we can measure S (of the macrostates).
we begin by stating that the total entropy is:
$$ S_{tot}= S + S_{micro}$$
where $S$ is the entropy that we can actually measure and $S_{micro}$ is the 
entropy associated with the microstates. It is the expected value of the variuos
microstate Entropy associated with each macrostate:
$$
S_{micro}=\langle S_{i}\rangle=\sum_i P_iS_i
$$
Where $S_i$ is given by Boltzmann formula $S_i=k_B\ln n_i$.  
Each macrostate has $n_i$ microstate and $\sum_i n_i=N $, the toatl number N of 
microstate of the system.
So the probability of finding the system in the i-th macrostate is:
$$
P_i= \frac {n_i} N
$$
Now we can derive the Gibbs formula, that express the entropy $S$ in terms of 
the probability $P_i$.  
$$
\begin{aligned}
S&= S_{tot} - S_{micro}\\
S&=k_B\left(\ln(N)-\sum_i P_i \ln(n_i)\right)\\
S&=k_B\left(\sum_i P_i(\ln(N) - \ln(n_i)) \right)\\
S&=-k_B\sum_i P_i \ln\left(\frac {n_i} N\right)\\
S&=-k_B\sum_i P_i \ln P_i
\end{aligned}
$$

# NOTE

### **Thermal equilibrium**

- Two or more systems are said to be in thermal equilibrium if their energy
  content and temperature are no longer changing with times.  
- Systems in thermal equilibrium have the same temperature
