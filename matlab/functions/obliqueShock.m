function beta = obliqueShock(theta,MachIn,n)
%OBLIQUESHOCK calc shock angle
%   INPUT:
%       theta: deviation angle [rad]
%       MachIn: Mach of the stream before the shock
%       n: degree of freedom of the gas, ex. for air n=5/2
%   OUTPUT:
%       beta: shock angle [rad]
%   NOTE:
%       initial guess of beta is equal to the limit shock angle arcsin(1/M) if
%       MachIn<4 or is given by small perturbation theory otherwise.
%       this assure convergence to the weak solution

    function y=func(beta,n)
        y=(1/(1+n))*(1+ n/(MachIn*MachIn*sin(beta)*sin(beta))) - tan(beta - theta)/tan(beta);
    end

if MachIn <=4 && MachIn > 1
    beta0 = asin(1/MachIn);
elseif MachIn > 4
    beta0 = theta*((n+1)/(2*n) + sqrt( ((n+1)/(2*n))^2 + (MachIn*theta)^-2 ));
end

beta = fzero(@(beta) func(beta,n), beta0);

end

