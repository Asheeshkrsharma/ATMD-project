%Author: Asheesh Sharma
%License: UIA (Use It Anywhere)
%Disclaimer: The UIA licesne applies to any part of this code except for
%any where it is specifically mentioned. You are given this software for
%free so dont try to sell it and do not bother me if some thing is broken.
%Finally, to all the students out there. The challenge is not to complete a
%task. It is to not plagiarize. :D 
%PS: This is based on a book called Orbital dynamics. You can find it for
%free online.
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
function [fdot, gdot] = fDot_and_gDot(x, r, ro, a,mu)
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
%
% This function calculates the time derivatives of the
% Lagrange f and g coefficients.
%
% mu - the gravitational parameter (kmˆ3/sˆ2)
% a - reciprocal of the semimajor axis (1/km)
% ro - the radial position at time t (km)
% t - the time elapsed since initial state vector (s)
% r - the radial position after time t (km)
% x - the universal anomaly after time t (kmˆ0.5)
% fDot - time derivative of the Lagrange f coefficient (1/s)
% gDot - time derivative of the Lagrange g coefficient
% (dimensionless)
%
% User M-functions required: stumpC, stumpS
% ------------------------------------------------------------
z = a*x^2;
%...Equation 3.66c:
fdot = sqrt(mu)/r/ro*(z*stumpS(z) - 1)*x;
%...Equation 3.66d:
gdot = 1 - x^2/r*stumpC(z);
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜