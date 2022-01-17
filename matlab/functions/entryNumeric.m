function [h,v,gamma,q,a] = entryNumeric(gas,ship,rho,T,g, h0, v0, gamma0)
% ENTRYNUMERIC use ode45 to integrate balistic entry equations
%   hp: gamma is not constant, B is constant
% INPUT:
%     ship: capsule struct
%         ship.B: balistic parameter (c_d*A/M)
%     rho: rho(h), function handler providing atmospheric model for the density 
%     g: g(h),  function handler providing gravity model
%     v0: initial velocity [m/s] 
%     gamma0: initial reentry angle [deg] 
% OUTPUT:
%     h: altitude vector [m]
%     v: velocity vector [m/s]
%     gamma: reentry angle vector [deg]
%     q: heat fluxes vector [Kw/m^2]
%     a: deceleration vector [g]
global r0

r0=6371e3; % [m]
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

[h,x]=ode45(@(h,x) balistic_rhs(h,x, ship, g, rho),[h0 0],[v0,gamma0],options);

v=x(:,1);
gamma=rad2deg(x(:,2));

%calc heat fluxes and deceleration
q= ship.St *rho(h).*v.*((v.^2)/2 + gas.cp*(T(h) - ship.Tw));

a = -0.5*ship.B*rho(h).*(v.^2) + transpose(g(h)).*sin(x(:,2));


end

