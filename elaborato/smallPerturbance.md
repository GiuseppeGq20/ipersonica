---
title: "small perturbance theory" 
author: "Giuseppe Giaquinto" 
---

# Mach indipendence

in the limit of hypersonic flow  and small perturbance the Mach number past a
shock wave does not depend on the upstream Mach number,infact: 

$$
M_2\rightarrow \frac 1 {n+1} \frac 1 { (\frac {n+1}n -1)^2} \frac 1 {\theta^2}\\
M_1 \rightarrow \infty
$$

# Equation of motion

## Nondimensionalising of the equation of motion in the limit of samll perturbance

For the reference quantity to non-dimensionalize the EOF we the jump in the
thermo-fluiddynamic quantities past the shock-wave and we can then ignore al
the therms of orde $\tau^2$ or grater while keeping the terms proportional to
any power of $K=M_{\infty}\tau$, that is the *hypersonic symilarity parameter*,
if we carry out the procedure for the euler equation we get:  
$$
\begin{aligned}
&\frac{\partial \rho}{\partial x} + 
v\frac{\partial \rho v}{\partial y} +
w\frac{\partial \rho w}{\partial z} =0\\
&\frac{\partial u}{\partial x}+
v\frac{\partial u}{\partial y}+
w\frac{\partial u}{\partial z}= -\frac 1 \rho \frac{\partial p}
{\partial x}\\
&\frac{\partial v}{\partial x}+
v\frac{\partial v}{\partial y}+
w\frac{\partial v}{\partial z}= -\frac 1 \rho \frac{\partial p}
{\partial y}\\
&\frac{\partial w}{\partial x}+
v\frac{\partial w}{\partial y}+
w\frac{\partial w}{\partial z}= -\frac 1 \rho \frac{\partial p}
{\partial z}\\
& \frac{\partial}{\partial x}\left(\frac p {\rho^\gamma}\right)+
v\frac{\partial}{\partial y}\left(\frac p {\rho^\gamma}\right)+
w\frac{\partial}{\partial z}\left(\frac p {\rho^\gamma}\right)=0
\end{aligned}
$$  
and their are paired with the following bondary conditions:

- downstream shock thermo fluid dynamic conditions:  
  $$
  \begin{aligned}
  u&= 1 - \frac{n+1}{n}\theta^2\\
  v&=\theta\\
  w&=\theta\\
  \rho&= n+1\\
  p&=\frac{(n+2)(n+1)}{n^2}M_\infty^2\theta^2\\
  \end{aligned}
  $$  
- slip flow on the body:  
  $$n_x + vn_y+wn_z=0$$  

> Note: we can se that in this set of partial differential equations the one
> relative to $x$ component of the momentum is indipendent from the others, and
> the $u$ does not appear in the boundary conditions, thus we can ignore it and
> solve the system without it and then calculate the $u$ component of the
> velocity from the energy equation.  

## Analisys of the flow field of an airfoil in hypersonic regime in the limits of the small perturbance theory

for a 2D airfoil moving at $Ma>>1$ in the limit of small perturbance theory we 
can express the equation of motion as:  
$$
\begin{aligned}  
&\frac{\partial \rho}{\partial x} +\frac{\partial \rho v}{\partial y}=0\\
&\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y} =- \frac{1}{\rho} 
\frac{\partial p}{\partial y}\\
&\frac{\partial}{\partial x} \frac p {\rho^\gamma} + v \frac{\partial}{\partial y}
 \frac p {\rho^\gamma} = 0\\
\end{aligned}
$$  
We are assuming uniform far field boundary condition, so that the flow is Homoenthalpic.  
Furthermore, since we are assuming that the disturb induced by the airfoil in the
flow field we neglect the curvature of front shock-wave,so that the flow is also 
Homoentropic. This allows us to replace the isentropic equation  with:  

$$
\frac{p}{\rho^\gamma}=Z=cost
$$  

and Z can be practically calculated from the condition past the shockwave, $Z= \frac{p_2}{\rho_2^\gamma}$.  

### conservative form
We can further elaborate this equation, thanks to the homoentropicity hypotesis to
express them in a conservative form. In fact, we can derive the following:  

$$
\begin{aligned}
&p\nu^\gamma=Z \qquad \nu= \frac 1 \rho\\ 
&\frac{\partial pv^\gamma}{\partial y}=0\\
&p\frac{\partial \nu}{\partial y}= - \frac \nu \gamma \frac{\partial p}{\partial y}
\end{aligned}
$$  

and  

$$
\begin{aligned}
& \nu \frac{\partial p}{\partial y} =  \frac{\partial p \nu}{\partial y} -
p\frac{\partial \nu}{\partial y}\\
& \nu \frac{\partial p}{\partial y} =  \frac{\partial p \nu}{\partial y} +
\frac \nu \gamma \frac{\partial p}{\partial y}\\
& \nu \frac{\partial p}{\partial y} = \frac \gamma {\gamma - 1}\frac{\partial p \nu}{\partial y}
\end{aligned}
$$

finally we can rewrite the equation of motion:

$$
\begin{aligned}
&\frac{\partial \rho}{\partial x} +\frac{\partial \rho v}{\partial y}=0 \\
&\frac{\partial v}{\partial x} + \frac {\partial}{\partial y}\left( \frac {v^2} 2 
+ \frac \gamma {\gamma - 1}\frac p \rho \right)=0 \\
\end{aligned}
$$ (1)

> the fact that we could express this equation in the conservative form allows as
> to more easily apply a finite volume scheme to integrate them.  

We can than obtain the $u$ component of the velocity via the Energy equation:  

$$
u + V_\infty = 
\sqrt{2\left(H_\infty- \frac \gamma {\gamma-1}Z\rho^{\gamma-1} \right) - v^2}
$$

### Boundary conditions

### numerical scheme

## Results

# Hayes analogy

In the previous section we elaborate on the 2D EOF, but we can recognise that 
eqs. (1) can be interpreted as the 1D  instationary euler equation where the role
of the time variable is taken by $x$.
