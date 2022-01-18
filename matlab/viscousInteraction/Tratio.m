function T_ratio= Tratio(M,gas,Tw)
T_ratio=0.42 + 0.032*(M^2) +0.58*(Tw/gas.T);
end

