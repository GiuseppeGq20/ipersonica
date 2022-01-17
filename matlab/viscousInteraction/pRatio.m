function p_ratio = pRatio(gas,M,theta,method)
%pRatio calc pressure ratio for given method
%   Input:
% gas: gas like struct
% M: freestream Mach number
% theta: profile inclination angle distribution
% method: ["shock_expansion","tangent_wedge"], method by which compute pressure ratio
%   Output:
% p_ratio: pressure ratio

switch method
    
    case "tangent_wedge"
        
        p_ratio= 1 + gas.gamma*(M*theta) + gas.gamma*((gas.gamma+1)/4)*(M*theta).^2;
        
    case "shock_expansion"
        
        beta=obliqueShock(theta(1),M,gas.n);
        p2_pinf=(1/(gas.n+1))*((gas.n+2) * (M*sin(beta))^2 -1 );
        theta=theta-theta(1);
        
        p_ratio=p2_pinf*(1+ M*theta/gas.n).^(gas.n*gas.gamma);
        
end

end

