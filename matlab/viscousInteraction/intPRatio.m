function Ip_ratio= intPRatio(p_ratio,x)
% calc cumulative integral in eq (3.2) pag 361

Ip_ratio=cumtrapz(x(1:end-1),p_ratio);

end

