% flate plate hypersonic boundary layer using viscous interaction theory
clc;clear;close all;

%flow data
gas.gamma=1.4;
gas.n=5;
gas.T=235; %[K]
gas.Taw=2143; %[K]
gas.v=1840; %[m/s]
gas.p=574; %[N/m^2] 
gas.rho=8.46e-3; %[kg/m^3]
gas.R=287.5;
gas.cp=gas.gamma*gas.R/(gas.gamma-1);
gas.T0=(gas.cp*gas.T + (gas.v*2)/2)/gas.cp; %[K]

M=5.97;
Re_inf=9.87e5;

%plate data
L=1; %[m]
x=linspace(0,L,80);
T_w=gas.Taw;

% boundary layer parameters
C_w=1; %(rho_W*mu_w)/(rho_e*mu_e) asumed parameter, it reproduces the results from the text, 

Re_x=(Re_inf/L).*x; % Rxe/Re_inf=x/L
Chi=(M^3)*(C_w./Re_x).^0.5; % Chi is singular at zero

%viscous interaction loop
[p_ratio_hot,delta_star_hot]=viscousInteraction(Chi,x,M,T_w,gas,"tangent_wedge"); %hot wall

[p_ratio1_hot,delta_star1_hot]=viscousInteraction(Chi,x,M,T_w,gas,"shock_expansion"); %hot wall

[p_ratio_cold,delta_star_cold]=viscousInteraction(Chi,x,M,gas.T,gas,"tangent_wedge"); %cold wall

p_ratio_a_hot=pRatioAnalytic(gas,Chi,T_w); % p_ratio analytic

%calc skin friction coefficient
cf=skinFriction(gas,p_ratio1_hot,Re_x,M,T_w); %hot wall
cf_cold=skinFriction(gas,p_ratio_cold,Re_x,M,gas.T); %cold wall

% calc drag coefficient
cd_hot=trapz(x,[cf(2),cf(2:end)])/L %hot wall
cd_cold=trapz(x,[cf_cold(2),cf_cold(2:end)])/L %cold wall 

%% plotting

%plot deltaStar distribution T_w=T_aw
plot(x,[delta_star_hot,NaN])
xlabel("x")
ylabel("\delta^*")
grid()

%plot pressure distribution
figure()
plot(x,[p_ratio_hot,NaN],"displayname","numerical method")
hold on
plot(x,p_ratio_a_hot,"-.","displayname","analytical method")
xlabel("x")
ylabel("p/p_\infty")
legend()
grid()

%plot skin friction distribution
figure()
plot(x,cf)
ylabel("c_f")
xlabel("x")
grid()

% plot viscous interaction parameter
figure()
plot(x,Chi)
xlabel("x")
ylabel("\chi")
grid()

%plot comparison deltaStar
labels={'T_w=T_{aw}','T_w=T_{\infty}'};
figure()
plot(x,[delta_star_hot,NaN],"-")
hold on
plot(x,[delta_star_cold,NaN],"-.")
xlabel("x")
ylabel("\delta^*")
legend(labels)
grid()

%plot comparison p ratio
figure()
plot(x,[p_ratio_hot,NaN],"-")
hold on 
plot(x,[p_ratio_cold,NaN],"-.")
xlabel("x")
ylabel("p/p_\infty")
legend(labels)
grid()

%plot comparison cf
figure()
plot(x,cf,"-")
hold on 
plot(x,cf_cold,"-.")
xlabel("x")
ylabel("c_f")
legend(labels)
grid()