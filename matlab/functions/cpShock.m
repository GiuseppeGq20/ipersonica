function cp = cpShock(Ma,n,beta)
%%calc cp past a shock wave
%input:
% Ma: Mach number before the shock
% n: thermodynamic  dof of the fluid
% beta: shock angle in radians
%output:
% cp: pressure coefficient
if nargin==2
    beta=0;
end

cp= (2*n/(n+1)) * (sin(beta).^2 - Ma^-2);
end

