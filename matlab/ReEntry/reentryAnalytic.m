%% rientro balistico
%script che implementa soluzione analitica semplificata e
%delle equazioni del rientro atmosferico

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
gas.T0= 288 ; %[K]
gas.R= 287 ; %[J/kgK]
gas.gamma=1.4;

%% definizione struct capsula
ship.B= c_d*A/M;
ship.St=St;

%% soluzione analitica
% angolo di rientro costante
% atmosfera isoterma

h=linspace(0,h0);


gamma=[4,5,7,8];
Legend={};
for i=1:length(gamma)
   
    Legend{i}=strcat("\gamma=",num2str(gamma(i)),"°");
end

figure(2)
for i=1:length(gamma)
    [v, Ma, q, a]=entryAnalytic(h,g0,gas,ship,gamma(i),v0);
    plot(q/1000,h)
    hold on
end
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
ylabel("H [Km]")
xlabel("Ma_{\infty}")
legend(Legend)
hold off

%%nota
% sottostima flussi termici e accelerazioni rispetto al testo
