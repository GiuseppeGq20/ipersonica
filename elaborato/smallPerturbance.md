---
title: "small perturbance airfoil"  
author: "Giuseppe Giaquinto"  
---

# Analisys the flow field of an airfoil in hypersonic regime in the limits of the small perturbance theory

## Equation of motion
for a 2D airfoil moving at $Ma>>1$ in the limit of small perturbance theory we 
can express the equation of motion as:
$$
\begin{align*}  
&\frac{\partial \rho}{\partial x} +\frac{\partial pv}{\partial y}=0\\
&\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y} =- \frac{1}{\rho} 
\frac{\partial p}{\partial y}\\
&\frac{\partial}{\partial x} \frac p {\rho^\gamma} + v \frac{\partial}{\partial y}
 \frac p {\rho^\gamma} = 0\\
\end{align*}
$$
We are assuming uniform far field boundary condition, so that the flow is Homoenthalpic.  
Furthermore, since we are assuming that the disturb induced by the airfoil in the
flow field we neglect the curvature of front shock-wave,so that the flow is also 
Homoentropic. This allows us to replace the isentropic equation  with:
$$
\frac{p}{\rho^\gamma}=Z=cost
$$
and Z can be practically calculated from the condition past the shockwave, 
$Z= \frac{p_2}{\rho_2^\gamma}$.  
### conservative form
We can further elaborate this equation, thanks to the homoentropicity hypotesis to
express them in a conservative form. In fact, we can derive the following:
$$
\begin{align*}
&p\nu^\gamma=Z \qquad \nu= \frac 1 \rho\\ 
&\frac{\partial pv^\gamma}{\partial y}=0\\
&p\frac{\partial \nu}{\partial y}= - \frac \nu \gamma \frac{\partial p}{\partial y}
\end{align*}
$$
and
$$
\begin{align*}
& \nu \frac{\partial p}{\partial y} =  \frac{\partial p \nu}{\partial y} -
p\frac{\partial \nu}{\partial y}\\
& \nu \frac{\partial p}{\partial y} =  \frac{\partial p \nu}{\partial y} +
\frac \nu \gamma \frac{\partial p}{\partial y}\\
& \nu \frac{\partial p}{\partial y} = \frac \gamma {\gamma - 1}\frac{\partial p \nu}{\partial y}
\end{align*}
$$
finally we can rewrite the equation of motion:
$$
\begin{align*}
&\frac{\partial \rho}{\partial x} +\frac{\partial pv}{\partial y}=0\\
&\frac{\partial v}{\partial x} + \frac {\partial}{\partial y}\left( \frac {v^2} 2 
+ \frac \gamma {\gamma - 1}\frac p \rho \right)=0\\
\end{align*}
$$
>the fact that we could express this equation in the conservative form allows as
to more easily apply a finite volume scheme to integrate them.  

We can than obtain the $u$ component of the velocity via the Energy equation:
$$
u + V_\infty = \sqrt{
    2\left(H_\infty- \frac \gamma {\gamma-1}Z\rho^{\gamma-1} \right) - v^2
}
$$
### Boundary conditions



### numerical scheme

## Results
