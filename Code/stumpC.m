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
function c = stumpC(z)
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
%
% This function evaluates the Stumpff function C(z) according
% to Equation 3.50.
%
% z - input argument
% c - value of C(z)
%
% User M-functions required: none
% ------------------------------------------------------------
if z > 0
c = (1 - cos(sqrt(z)))/z;
elseif z < 0
c = (cosh(sqrt(-z)) - 1)/(-z);
else
c = 1/2;
end
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜