function [ E ] = true2eccentric( nu,ecen )
%Converts true anamoly to eccentric anomaly
E = 2*atan(sqrt((1-ecen)/(1+ecen))*tan(nu/2));
end