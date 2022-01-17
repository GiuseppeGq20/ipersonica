function [p_ratio,delta_star] = viscousInteraction(Chi,x,M,T_w,gas,method)

% viscous interaction loop
err=1;
theta=zeros(1, length(x)-1);
delta_star=zeros(1,length(x)-1);
dx=x(2:end)-x(1:end-1);
p_ratio=ones(1,length(x)-1);

%iter=0;
while err>1e-3

    temp=deltaStar(gas,p_ratio,Chi,M,x,T_w);
    
    theta=atan(diff(temp)./diff(x(1:end-1)));
    theta(1)=theta(2); % to eliminate singularities in p_ratio calculation
    theta=[theta,theta(end)];
    
    %tangent wedge
    switch method
        case "shock_expansion"
             p_ratio=pRatio(gas,M,theta,"shock_expansion");
        case "tangent_wedge"
            p_ratio=pRatio(gas,M,theta,"tangent_wedge");
    end
    %add shock expansion method
    
    err=norm((delta_star-temp));
    delta_star=temp;
    
    %iter=iter+1;
end
p_ratio(1)=NaN;
end

