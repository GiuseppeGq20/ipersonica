function cf= skinFrictionm(gas,p_ratio,Re_x,M,T_w)
 %%SKINFRICTION calc skin friction distribution on a flat plate
 % Input:
% gas: gas like struct
% p_ratio: pressure ratio
% Re_x: local Reynolds array
% M: free stream Mach number
% T_w: wall temperature
%Output:
%cf: skin friction coefficient

% calc skin friction coefficient you, must use Ve
% suppose Te=T_inf=cost

Ve=(0.5*(gas.v^2 + gas.R*gas.T*log(p_ratio))).^.5;
Ve(1)=Ve(2);
Re_x=Re_x.*([Ve,Ve(end)]/gas.v);
cf=0.664*(Tratio(M,gas,T_w)^(-1/3))./(Re_x.^0.5);
 
end
