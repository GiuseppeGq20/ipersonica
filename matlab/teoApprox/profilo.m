%% Esercitazione teorie approssimate profilo

clear all; clc; close all;
%condizioni moto
Ma=8;
alpha= deg2rad(8); % angolo d' attacco
n=5; % dof gas
gamma=1.4;
%% geometria
%profilo starfighter
c=1;
tau=0.15;
R= (1/tau + tau)*c/4; % radius of curvature is constant

angle=linspace(pi,0,150);
x=(1+cos(angle))/2;
yDorso=-(R- 0.5*tau) + (R^2 - (x- 0.5).^2).^0.5;
yVentre= +(R- 0.5*tau) - (R^2 - (x- 0.5).^2).^0.5;


deltayDorso= yDorso(2:end) - yDorso(1:end-1);
deltayVentre= yVentre(2:end) - yVentre(1:end-1);

deltax=x(2:end)-x(1:end-1);
detalDorso = realsqrt(deltax.^2 + deltayDorso.^2); detalDorso=[detalDorso,NaN];
deltalVentre = realsqrt(deltax.^2 + deltayVentre.^2); deltalVentre=[deltalVentre,NaN];

thetaDorsoGeom=atan2(deltayDorso,deltax); 
thetaDorsoGeom=[thetaDorsoGeom,NaN];

thetaDorso=thetaDorsoGeom - alpha;

thetaVentreGeom=atan2(deltayVentre,deltax);
thetaVentreGeom=-[thetaVentreGeom,NaN]; % qui c'è un segno meno per avere segni consistenti dopo

thetaVentre= thetaVentreGeom + alpha;

%% cp Newton

cpDorsoN=2*sin( thetaDorso(thetaDorso>0)).^2; cpDorsoN=[cpDorsoN,nan(1,length(x)-length(cpDorsoN))]; % orlo per plottare
cpVentreN=2*sin( thetaVentre(thetaVentre>0) ).^2; cpVentreN=[cpVentreN,nan(1,length(x)-length(cpVentreN))]; % orlo per plottare

% calcolo cl e cd
[ClN,CdN]=CalcClCd(alpha,thetaDorsoGeom,thetaVentreGeom,detalDorso,deltalVentre,cpDorsoN,cpVentreN);

%% cp Busemann (nell' approssimazione che la curvatura sia molto piccola)
cpDorsoB=cpDorsoN - yDorso*2/R; cpDorsoB(cpDorsoB<0)= NaN; % cp negativi non hanno senso nella teoria di Newton
cpVentreB=cpVentreN + yVentre*2/R; cpVentreB(cpVentreB<0) = NaN;

% calcolo Cl e Cd
[ClB,CdB]=CalcClCd(alpha,thetaDorsoGeom,thetaVentreGeom,detalDorso,deltalVentre,cpDorsoB,cpVentreB);

%% cp cono tangente
cpDorsoC= cpCuneoTangente(Ma,n,thetaDorso);
%correction for curvature
%  cpDorsoC= cpDorsoC - (2*yDorso)/R;

cpVentreC = cpCuneoTangente(Ma,n,thetaVentre);
% cpVentreC= cpVentreC + (2*yVentre)/R;

%% cp urto espansione
cpDorsoUE=cpUrtoEspansione(Ma,thetaDorso,n);

cpVentreUE=cpUrtoEspansione(Ma,thetaVentre,n);

%% plotting

% cpDorso
figure(1)

plot(x,cpDorsoN,'k-', 'DisplayName','Newton')
hold on
plot(x,cpDorsoB,'b-','DisplayName','Buseman')
hold on
plot(x,cpDorsoC,'r-.','DisplayName','Cono Tangente')
hold on
plot(x,cpDorsoUE,'g-.','linewidth',1,'DisplayName','Urto-espansione')
hold off
title("C_p sul dorso del profilo")
xlabel("x") 
ylabel("C_p")
legend
grid on

%cpVentre
figure(2)
plot(x,cpVentreN,'k-', 'DisplayName','Newton')
hold on
plot(x,cpVentreB,'b-','DisplayName','Buseman')
hold on
plot(x,cpVentreC,'r-.','DisplayName','Cono Tangente')
hold on
plot(x,cpVentreUE,'g-.','linewidth',1,'DisplayName','Urto-espansione')
hold off
title("C_p sul ventre del profilo")
xlabel("x")
ylabel("C_p")
legend
grid on


