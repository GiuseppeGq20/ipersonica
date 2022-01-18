%% Analytical and Numerical entry comparison
clc;clear ; close all;
%% dati capsula
M= 2800; %[kg]
d=2.8; %[m]
A=pi * (d/2)^2 ; %[m^2]
c_d=1;
St=0.01;

%% dati iniziali
h0=120e3; %[m]
v0= 7500; %[m/s]
g0=9.81; %[m/s^2]

%% dati gas atm
gas.rho0= 1.225; %[kg/m^3]
gas.T0= 231.5 ; %[K]
gas.R= 287 ; %[J/kgK]
gas.gamma=1.4;
gas.cp=1004;
%% definizione struct capsula
ship.B= c_d*A/M;
ship.St=St;
ship.Tw=900;

%% dati modello di atmosfera
atm=importdata("USSA76.dat");
atmH=atm.data(:,1)*1000;
atmRho=atm.data(:,4);
atmT=atm.data(:,3);
%% comparison

%numeric solution
gamma0=4; %[deg]
r0=6371e3; % [m]

rho=@(x) interp1(atmH,atmRho,x);
g=@(x) g0*(1/(1 + (x/r0)).^2);
T=@(x) interp1(atmH,atmT,x);
[h_n,v_n,gamma_n,q_n,a_n] = entryNumeric(gas,ship,rho,T,g, h0, v0, deg2rad(gamma0));


h=linspace(0,h0);
gamma=[4,5,7,8];
Legend={};
for i=1:length(gamma)
   
    Legend{i}=strcat("\gamma=",num2str(gamma(i)),"°");
end
Legend{length(gamma)+1}='numeric 4°';

%plotting
figure(2)
for i=1:length(gamma)
    [v, Ma, q, a]=entryAnalytic(h,g0,gas,ship,gamma(i),v0);
    plot(q/1000,h)
    hold on
end
plot(q_n/1000,h_n,'k--')
ylabel("H [Km]")
xlabel("q [Kw/m^2]")
legend(Legend)
hold off

figure(3)
for i=1:length(gamma)
    [v, Ma, q, a]=entryAnalytic(h,g0,gas,ship,gamma(i),v0);
    plot(v,h)
    hold on
end
plot(v_n,h_n,'k--')
ylabel("H [Km]")
xlabel("v [m/s]")
legend(Legend)
hold off

figure(4)
for i=1:length(gamma)
    [v, Ma, q, a]=entryAnalytic(h,g0,gas,ship,gamma(i),v0);
    plot(-a,h)
    hold on
end
plot(-a_n/g0,h_n,'k--')
ylabel("H [Km]")
xlabel("a/g_0")
legend(Legend)
hold off

figure(5)
for i=1:length(gamma)
    [v, Ma, q, a]=entryAnalytic(h,g0,gas,ship,gamma(i),v0);
    plot(Ma,h)
    hold on
end
plot(v_n/sqrt(gas.gamma*gas.R*gas.T0),h_n,'k--')
ylabel("H [Km]")
xlabel("Ma_{\infty}")
legend(Legend)
hold off
