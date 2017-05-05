%Author: Asheesh Sharma
%License: UIA (Use It Anywhere)
%Disclaimer: The UIA licesne applies to any part of this code except for
%any where it is specifically mentioned. You are given this software for
%free so dont try to sell it and do not bother me if some thing is broken.
%Finally, to all the students out there. The challenge is not to complete a
%task. It is to not plagiarize. :D 
%PS: This is based on a book called Orbital dynamics. You can find it for
%free online.
function [ a,e,incl,Omega,w,nu ] = RV2COE( Rijk, Vijk, mu )
%{
This function takes the state vector (Rijk and Vijk)
and computes the classical orbital elements

Inputs:
    Rijk - Position vector in Earth Centered Inertial (km)
    Vijk - Velocity vector in Earth Centered Inertial (km)

Outputs:
    a     - semi-major axis of the orbit (km)
    e     - eccentricity
    i     - inclination (rad)
    Omega - right ascension of the ascending node (rad)
    w     - argument of perigee (rad)
    nu    - true anomaly (rad)
%}
eps = 1.e-10;
%mu=398600.4418; %Earth

r = norm(Rijk);
v = norm(Vijk);
vr = dot(Rijk,Vijk)/r;
H = cross(Rijk,Vijk);
h = norm(H);
%Get inclination
incl = acos(H(3)/h);
%Calculate normal to angular velocity
N = cross([0 0 1],H);
n = norm(N);
%Calculate right ascension of the ascending node
if n ~= 0
    Omega = acos(N(1)/n);
    if N(2) < 0
        Omega = 2*pi - Omega;
    end
else
    Omega = 0;
end
%Find eccentricity vector and get eccentricity
E = 1/mu*((v^2 - mu/r)*Rijk - r*vr*Vijk);
e = norm(E);
%Calculate argument of perigee
if n ~= 0
    if e > eps
        w = acos(dot(N,E)/n/e);
        if E(3) < 0
            w = 2*pi - w;
        end
    else
        w = 0;
    end
else
    w = 0;
end
%Calculate true anomaly
if e > eps
    nu = acos(dot(E,Rijk)/e/r);
    if vr < 0
        nu = 2*pi - nu;
    end
else
    cp = cross(N,Rijk);
    if cp(3) >= 0
        nu = acos(dot(N,Rijk)/n/r);
    else
        nu = 2*pi - acos(dot(N,Rijk)/n/r);
    end
end
%.Calculate semimajor axis
a = h^2/mu/(1 - e^2);
end

