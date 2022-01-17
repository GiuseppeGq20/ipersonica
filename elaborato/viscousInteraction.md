---
title: "Viscous Interaction"
author: "Giuseppe Giaquinto"
---

# Viscous interaction

One possible way to compute boundary layers parameters is the viscous interaction
method.  
It is an approximate method, because it relays on integral parameters ($\delta^*$),
and to calculate the skin friction coefficient we must use empirical formulas.  
The name viscous interaction refers to the fact that by the apparent thickness $\delta^*$ 
the boundary layer interact with pressure field outside the boundary layer, and thus,
it may or may not affect it enough to influence the boundary layer itself.  
The parameter that give us a mesure of the strength of the viscous interaction is
$\Chi$, defined as:
$$
\Chi = M_\infty^3 \sqrt{\frac{C_w}{Re_\infty}}
$$

Where $C_w=\rho_w \mu_w/\rho_e \mu_e \simeq 1$.  
We can have three cases:

- $\Chi < 1$, there is no viscous interaction, the pressure field outside the 
boundary layer is independet from the boundary layer itself
- $\Chi \simeq 1$, weak viscous interaction, the pressure field can still be calculated independetly from the boundary layer, but we must recalculate it after
the "apparent geometry" is known
- $\Chi > 1$, strong viscous interaction, the boundary layer and the pressure field
must be resolved simultaneously

## numerical method

We can define a numerical procedure, that can handle bot weak and strong interactions.
It is summarized in these steps:

1) compute the pressure field in the non viscous flow
2) compute $\delta^*$ distribution on the body
3) update body geometry with the $\delta^*$ distribution
4) repeat steps 1), 2), 3) until convergence

## case study: flate plate

To solve the boundary layer of an hypersonic flat plate with viscous interaction
with the aid of the ...

### A formula for the $\delta^*$ calculation

A formula for the calculation of $\delta^*$ obtained by Cox and Crabtree 
(suitable for use in a viscous interaction process) is given here without proof:
$$
\frac{\delta^*}{x}= \frac{\gamma -1}{\gamma +1} (0.664 + 1.73 \frac{T_w}{T_0})
\frac{\Chi}{M_\infty}\frac{p_\infty}{p}\frac{1}{\sqrt{x}}\left( \int_0^x \frac{p}{p_\infty} dx \right)^{\frac{1}{2}}
$$

To get an estimate of the skin friction coefficient $c_f$ we can use the blasius
formula, but we evaluate the thermodynamic properties at $T=T^*$, where $T^*$ 
is the reference temperature.
$$
C_f=\frac{0.664}{\sqrt{Re_{xe}}}
$$

### Results