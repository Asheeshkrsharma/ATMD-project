function  [Rijk, Vijk] = COE2RV(a, e, incl, Omega, w, nu, mu)
%{
This function computes takes the classical orbital elements
and computes the state vector

Inputs:
    a     - semi-major axis of the orbit (km)
    e     - eccentricity
    i     - inclination (rad)
    Omega - right ascension of the ascending node (rad)
    w     - argument of perigee (rad)
    nu    - true anomaly (rad)

Outputs:
    Rijk - Position vector in Earth Centered Inertial (km)
    Vijk - Velocity vector in Earth Centered Inertial (km)
%}
%mu=398600.4418; %Earth
h=sqrt(mu*a*(1-e^2));

%Calculate the position and velocity vectors in perifocal frames
rp = (h^2/mu) * (1/(1 + e*cos(nu))) * (cos(nu)*[1;0;0] + sin(nu)*[0;1;0]); 
vp = (mu/h) * (-sin(nu)*[1;0;0] + (e + cos(nu))*[0;1;0]); 
  
%Create rotation matricies for 3-1-3 rotation
R3_W = [ cos(Omega)  sin(Omega)  0 
        -sin(Omega)  cos(Omega)  0 
            0        0     1]; 
  
R1_i = [1       0          0 
        0   cos(incl)  sin(incl) 
        0  -sin(incl)  cos(incl)]; 
  
R3_w = [ cos(w)  sin(w)  0  
        -sin(w)  cos(w)  0 
           0       0     1]; 
  
%Create full rotation matrix
Q_pX = (R3_w*R1_i*R3_W)'; 
%Calculate instantaneous postiton and velocity vectors
Rijk = Q_pX*rp; 
Vijk = Q_pX*vp; 
end

