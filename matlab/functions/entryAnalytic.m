function [v, Ma, q, a]=entryAnalytic(h,g0,gas,ship,gamma,v0)
% calcola le grandezze di rientro sfruttando il modello analitico:
% note:
%   atmosfera isoterma
%   angolo di rientro costante
%
% Input:
%     h: dimensional height vector in Km
%     g: gravity acceleration
%     gas: struct
%         gas.R: gas constant
%         gas.rho: initial gas density
%         gas.T0: stagnation gas temperature
%     ship: struct
%         ship.B: balistic parameter
%         ship.St: stanton mean number of the flow
%     gamma: initial entry angle in [deg]
%     v0: initial ship velocity
%
% Output:
%   v: velocity
%   Ma: mach number of the capsule
%   q: heat flux incoming in the capsule
%   a: capsule deceleration [g]

    
    beta=g0/(gas.R*gas.T0); % c'è un fattore 10^3 perchè ho espresso h in km
    
    z=beta*h;
    
    esp=exp(-z);
    
    alpha=(0.5*gas.rho0*ship.B)/(beta*sind(gamma));
    
    v=v0*exp(-alpha*esp);
    
    q=0.5 * gas.rho0 * ship.St * esp.*(v.^3);
    
    a= -alpha*beta*sind(gamma)* (1/g0) * (v.^2) .* esp;
    
    Ma=v./sqrt(gas.gamma*gas.R*gas.T0);
end
