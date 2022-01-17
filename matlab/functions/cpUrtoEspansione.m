function cp=cpUrtoEspansione(Ma,theta,n)
%% calc cp with shcok expansion theory
%INPUT:
% Ma: upstream mach
% n: gas d.o.f.
% theta: current angle deviation array
%OUTPUT:
%cp: cp distribution


%calc leading edge shock angle
beta= obliqueShock(theta(1),Ma,n); 

%calc downstream  normal Mach
M2=sqrt(...
    ((Ma*sin(beta))^2 + n) / ((n+2)* ((Ma*sin(beta))^2) -1) ...
    );
%calc shock pressure ratio
P=(1/(n+1))*((n+2)*(Ma*sin(beta))^2 - 1);

%calc downstream mach
M2=M2/sin(beta-theta(1));

% calc subsequent current deviation
theta=theta-theta(1);

%calc whole pressure ratio
PPinf=P * (1 + ((M2*theta)/n)).^(n+2);

% PPinf=zeros(1,length(theta));
% PPinf(1)=P;
% for i=2:length(PPinf)
%     
%     [M,nu,mu]=flowprandtlmeyer(gamma,M2);
%     nu-rad2deg(theta(i));
%     M2=flowprandtlmeyer(gamma,nu+rad2deg(-theta(i)),'nu');
%     PPinf(i)=PPinf(i-1)*(1 +(M2*theta(i))/n)^(n*gamma);
% 
%     
% end

%calc pressure coefficient
cp=(PPinf -1)*(2*n)/((n+2) * Ma^2);

end