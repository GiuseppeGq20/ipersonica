%% rientro balistico
%%script che implementa la soluzione
%%numerica delle equazioni del rientro atmosferico

clc;clear all; close all;

global r0

%% dati capsula
M= 2800; %[kg]
d=2.8; %[m]
A=pi * (d/2)^2 ; %[m^2]
c_d=1;
St=0.01;

%% dati iniziali
h0=120e3; %[km]
v0= 7500; %[km/s] ?
gamma0=deg2rad(4) ; %[deg]
g0=9.81; %[m/s^2]
r0=6371e3; % [m]
%% dati gas
gas.rho0= 1.225; %[kg/m^3]
gas.T0= 288 ; %[K]
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
%% soluzione numerica
% rho=@(h)exp(-h*g0/(gas.R*gas.T0));
rho=@(h) interp1(atmH,atmRho,h);
g=@(h)g0*(1/(1 + (h/r0)).^2);
T=@(h)interp1(atmH,atmT,h);

[h,v,gamma,q,a] = entryNumeric(gas,ship,rho,T,g, h0, v0, gamma0);

plot(gamma,h)
