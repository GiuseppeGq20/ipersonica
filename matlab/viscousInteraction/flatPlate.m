% flate plate hypersonic boundary layer using viscous interaction theory

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
x=linspace(0,L,40);
T_w=gas.Taw;

% boundary layer parameters
C_w=1; %(rho_W*mu_w)/(rho_e*mu_e) asumed parameter, it reproduces the results from the text, 

Re_x=(Re_inf/L).*x; % Rxe/Re_inf=x/L
Chi=(M^3)*(C_w./Rex).^0.5; % Chi is singular at zero

% viscous interaction loop
err=1;
theta=zeros(1, length(x)-1);
delta_star=zeros(1,length(x)-1);
dx=x(2:end)-x(1:end-1);
p_ratio=ones(1,length(x)-1);

iter=0;
while err>1e-3

    temp=deltaStar(gas,p_ratio,Chi,M,x,T_w);
    
    theta=atan(diff(temp)./diff(x(1:end-1)));
    theta=[theta,theta(end)];
    
    %tangent wedge
    p_ratio= 1 + gas.gamma*(M*theta) + gas.gamma*((gas.gamma+1)/4)*(M*theta).^2;
    
    %add shock expansion method
    
    err=norm((delta_star-temp));
    delta_star=temp;
    
    iter=iter+1
end

% calc skin friction coefficient you must use Ve
cf=0.664*(Tratio(M,gas,T_w)^(-1/3))./(Re_x.^0.5); %adjust cf calulation, see pag, 344

% calc drag coefficient
cd=trapz(x(2:end),cf(2:end))/L

%% plotting
plot(x,[delta_star,NaN])
figure()
plot(x,[p_ratio,NaN])
figure()
plot(x,cf)
figure()
plot(x,Chi)