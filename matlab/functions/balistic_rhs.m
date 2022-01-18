function xdot= balistic_rhs(h,x, ship, g,rho)
        
        global r0
        %x1 is v
        %x2 is gamma
        xdot=zeros(2,1);
        
        xdot(1) = 0.5* rho(h)*ship.B*x(1)/sin(x(2)) - g(h)/x(1);
        xdot(2) = 1/((r0 + h)*sin(x(2)))  - g(h)/((x(1)^2)*tan(x(2)));

end