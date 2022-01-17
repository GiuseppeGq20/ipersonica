% flate plate hypersonic boundary layer using viscous interaction theory
clc;clear;close all
%flow data
gas.gamma=1.4;
gas.n=5;
gas.T=235;
gas.Taw=2143;
gas.v=1840;
gas.p=574;
gas.rho=8.46e-3;
gas.R=287.5;
gas.cp=gas.gamma*gas.R/(gas.gamma-1);

gas.T0=(gas.cp*gas.T + (gas.v*2)/2)/gas.cp;

M=5.97;
Re_inf=9.87e5;

%plate data
L=1;
x=linspace(0,L,80);
T_w=gas.Taw;

% boundary layer parameters
C_w=1; %(rho_W*mu_w)/(rho_e*mu_e) asumed parameter, it reproduces the results from the text, 

Re_x=(Re_inf/L).*x; % Rxe/Re_inf=x/L
Chi=(M^3)*(C_w./Re_x).^0.5; % Chi is singular at zero

%viscous interaction loop
[p_ratio,delta_star]=viscousInteraction(Chi,x,M,T_w,gas,"tangent_wedge");

[p_ratio1,delta_star1]=viscousInteraction(Chi,x,M,T_w,gas,"shock_expansion");

% calc skin friction coefficient you, must use Ve
% suppose Te=T_inf=cost
Ve=(0.5*(gas.v^2 + gas.R*gas.T*log(p_ratio1))).^.5;
Ve(1)=Ve(2);
Re_x=Re_x.*([Ve,Ve(end)]/gas.v);

cf=0.664*(Tratio(M,gas,T_w)^(-1/3))./(Re_x.^0.5); %adjust cf calulation, see pag, 344

% calc drag coefficient
cd=trapz(x,[cf(2),cf(2:end)])/L

%% plotting

plot(x,[delta_star,NaN])

figure()
plot(x,[p_ratio,NaN])
figure()
plot(x,cf)
ylabel("c_f")
figure()
plot(x,Chi)