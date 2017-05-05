%Author: Asheesh Sharma
%License: UIA (Use It Anywhere)
%Disclaimer: The UIA licesne applies to any part of this code except for
%any where it is specifically mentioned. You are given this software for
%free so dont try to sell it and do not bother me if some thing is broken.
%Finally, to all the students out there. The challenge is not to complete a
%task. It is to not plagiarize. :D 
%PS: This is based on a book called Orbital dynamics. You can find it for
%free online.
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
