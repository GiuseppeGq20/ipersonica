function p_ratio= pRatioAnalytic(gas,Chi,Tw)
 
 p_ratio= 1 + (gas.gamma/2)*((gas.gamma-1)/(gas.gamma +1))*(0.664 + 1.73*(Tw/gas.T0)).*Chi;
 
end
