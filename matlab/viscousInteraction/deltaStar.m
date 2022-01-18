function delta = deltaStar(gas,p_ratio,Chi,M,x,Tw)
%DeltaStar calc delta star distribution along the plate
%   Detailed explanation goes here

%to assure no array dimension errors
Chi=Chi(1:end-1);

delta=((gas.gamma-1)/(gas.gamma+1))*(0.664 + 1.73*(Tw/gas.T0))*...
    (Chi/M).* p_ratio .* (intPRatio(p_ratio,x).* x(1:end-1)).^0.5;

end

