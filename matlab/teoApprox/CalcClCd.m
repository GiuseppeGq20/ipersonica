function [Cl,Cd]=CalcClCd(alpha,thetaDorso,thetaVentre,lDorso,lVentre,cpDorso,cpVentre)
%% calcola Cl e cd profilo 2D
%input:
% alpha: angolo d'attacco
% thetaDorso: angoli geometrici pannelli dorso
% thetaVentre: angoli geometrici pannelli ventre
% cpDorso: distribuzione cp dorso
% cpVentre: distribuzione cp ventre
%output:
% [Cl,Cd]: coefficienti di portanza e resistenza

Cx= - cpDorso.*lDorso.*cos(thetaDorso + pi/2) - cpVentre.*lVentre.*cos(thetaVentre + pi/2);
Cx(Cx==NaN)=0;
%Cx=sum(Cx,'omitnan');
Cx=sum(Cx);

Cy= - cpDorso.*lDorso.*sin(thetaDorso + pi/2) + cpVentre.*lVentre.*sin(thetaVentre + pi/2);
Cy(Cy==NaN)=0;
%Cy=sum(Cy,'omitnan');
Cy=sum(Cy);

Cd= Cx*cos(alpha) + Cy*sin(alpha);
Cl=-Cx*sin(alpha) + Cy*cos(alpha);
end