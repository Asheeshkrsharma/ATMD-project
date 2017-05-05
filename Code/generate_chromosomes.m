
%Author: Asheesh Sharma
%License: UIA (Use It Anywhere)
%Disclaimer: The UIA licesne applies to any part of this code except for
%any where it is specifically mentioned. You are given this software for
%free so dont try to sell it and do not bother me if some thing is broken.
%Finally, to all the students out there. The challenge is not to complete a
%task. It is to not plagiarize. :D 
function chrome = generate_chromosomes(num)
%Input num: Number of chromosomes to generate
%**********************Mission Time frame (Since epoch)********************
dt0 = (1000-108).*rand(num,1) + 108; %Time frame for completing transit orbit 1; leo1->EROS [days] Min 108 max 1498
dt1 = (1000-92).*rand(num,1) + 92; %Time frame for completing transit orbit 2; EROS->leo2 [days] Min 92 max 1500
a_c_leo1=(6.763353215358037E+03-2.063997893546698E+02).*rand(num,1) + 2.063997893546698E+02; %Semi-major axis, a (km) for circular orbits, it is the radius
nu_c_leo1=rad2deg((2*pi).*rand(num,1)); %True anomaly, nu (rad)
a_c_leo2=(6.763353215358037E+03-2.063997893546698E+02).*rand(num,1) + 2.063997893546698E+02; %Semi-major axis, a (km) for circular orbits, it is the radius
nu_c_leo2=rad2deg((pi).*rand(num,1)); %True anomaly, nu (rad)
incl_c_leo1=(pi.*rand(num,1)); %Inclination w.r.t XY-plane, i (rad)  
incl_c_leo2=(pi.*rand(num,1)); %Inclination w.r.t XY-plane, i (rad)  
chrome=[dt0 dt1 a_c_leo1 nu_c_leo1 a_c_leo2 nu_c_leo2 incl_c_leo1 incl_c_leo2];
end
