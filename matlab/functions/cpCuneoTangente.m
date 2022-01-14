function cp = cpCuneoTangente(Ma,n,theta)
%% calc cp distribution using tangent cone theory
%INPUT:
% Ma: upstream mach
% n: gas d.o.f.
% theta: current angle deviation vector

cp=NaN(1,length(theta));

    for i=1:length(theta)
    
        if theta(i)>=0
            beta=obliqueShock(theta(i),Ma,n);
            cp(i)= ((2*n)/(n+1))*((sin(beta)^2) - Ma^-2);
        else
            p2p=(1+ (Ma*theta(i))/n)^(n+2);
            cp(i)=(p2p-1)*(2*n)/((n+2)* Ma^2);
        end
    end
    
end