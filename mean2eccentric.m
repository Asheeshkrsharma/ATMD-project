function E = mean2eccentric(Me,ecen)
%Converts mean anamoly to eccentric anamoly
if Me < pi
    E = Me + ecen/2;
else
    E = Me - ecen/2;
end
En=1e10;
while abs(En-E) > 1e-8
    En = E + (Me - E + ecen*sin(E))/(1-ecen*cos(E));
    E=En;
end
end
