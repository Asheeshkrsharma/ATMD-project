%Author: Asheesh Sharma
%License: UIA (Use It Anywhere)
%Disclaimer: The UIA licesne applies to any part of this code except for
%any where it is specifically mentioned. You are given this software for
%free so dont try to sell it and do not bother me if some thing is broken.
%Finally, to all the students out there. The challenge is not to complete a
%task. It is to not plagiarize. :D 
%PS: This is based on a book called Orbital dynamics. You can find it for
%free online.
function [total,C3_total,vector]=evaluate_chromosome(chromosome)
dt0=chromosome(1);
dt1=chromosome(2);
a_c_leo1=chromosome(3);
nu_c_leo1=chromosome(4);
a_c_leo2=chromosome(5);
nu_c_leo2=chromosome(6);
incl_c_leo1=chromosome(7);
incl_c_leo2=chromosome(8);
%% Set parameters
%***********************************constants******************************
mu_sun=132712440018; % Keplerian Gravitaitonal parameter for bodies orbiting Sun [km^3/s^2]
mu_e=3.986004415e5; % Keplerian Gravitaitonal parameter for bodies orbiting Earth [km^3/s^2]

Re=6378.1363; % Radius of Earth [km]
%**************************************************************************


%**********************Orbital Elements************************************
%     Data queried from 
%     Solar System Dynamics Group, Horizons On-Line Ephemeris System
%     4800 Oak Grove Drive, Jet Propulsion Laboratory
%     Pasadena, CA  91109   USA
%     Information: http://ssd.jpl.nasa.gov/
%     Connect    : telnet://ssd.jpl.nasa.gov:6775  (via browser)
%                  http://ssd.jpl.nasa.gov/?horizons
%                  telnet ssd.jpl.nasa.gov 6775    (via command-line)
%-----------------------------------------------------------------------------
% Ephemeris / WWW_USER Thu Apr 27 08:05:11 2017 Pasadena, USA      / Horizons    
%-----------------------------------------------------------------------------
% Target body name: Earth (399)                     {source: DE431mx}
% Center body name: Sun (10)                        {source: DE431mx}
% Center-site name: BODY CENTER
r_earth=149598023; %Heliocentric Semi-major axis of earth [km]
rp_earth=147098290; %Heliocentric Perihelion of earth [km]
ra_earth=152098232; % Heliocentric Aphelion of earth [km]
a_e=(rp_earth+ra_earth)/2; %Average orbital radius of earth [km]
ecen_e=1.629513433577725E-02; %Eccentricity, e 
OM_e=348.73936*pi/180; %Longitude of Ascending Node, OM(rad)  
w_e=114.20783*pi/180; %Argument of Perifocus, w (rad)
me=5.972e24; %Mass of Earth [kg] (Gravitational perturbation not included)
nu_e=1.108533117592837E+02*pi/180; %True anomaly, nu (rad)
n_e=sqrt(mu_sun/a_e^3); %mean motion rad/s
incl_e=2.961070760988626E-03*pi/180; %Inclination w.r.t XY-plane,(rad)

%------------------------------Orbital Elements Eros----------------------------
% Ephemeris / WWW_USER Thu Apr 27 08:07:44 2017 Pasadena, USA      / Horizons    
%------------------------------Orbital Elements Eros----------------------------
% Target body name: 433 Eros (1898 DQ)              {source: JPL#610}
% Center body name: Sun (10)                        {source: DE431}
% Center-site name: BODY CENTER
%------------------------------Orbital Elements Eros----------------------------
a_d=2.180969722068203E+08; %Semi-major axis, a (km) 
ecen_d=2.225443543150867E-01; %Eccentricity, e 
incl_d=1.082785072634908E+01*pi/180; %Inclination w.r.t XY-plane, i (rad)  
w_d=1.788116740186882E+02*pi/180; %Argument of Perifocus, w (rad)
OM_d=3.043259456708146E+02*pi/180;  %Longitude of Ascending Node, OM(rad)  
nu_d=3.575782637883727E+02*pi/180; %True anomaly, nu (rad)
n_d=sqrt(mu_sun/a_d^3); %mean motion rad/s
%**************************************************************************

%------------------------------Orbital Elements for LEO1 Craft---------------------------
%*******************************************************************************
% Ephemeris / WWW_USER Thu Apr 27 13:13:53 2017 Pasadena, USA      / Horizons    
%*******************************************************************************
% Target body name: Craft (NA)                      {source: DE431mx}
% Center body name: Earth (399)                     {source: DE431mx}
% Center-site name: BODY CENTER
%*******************************************************************************
ecend_c_leo1=0; %Eccentricity, e for circular orbits is zero
OM_c_leo1=0;%Longitude of Ascending Node, OM(rad)    
w_c_leo1= 0; %Argument of Perifocus, w (rad)
%-------------------------Orbital Elements Craft---------------------------

%------------------------------Orbital Elements for LEO2 Craft---------------------------
%*******************************************************************************
% Ephemeris / WWW_USER Thu Apr 27 13:13:53 2017 Pasadena, USA      / Horizons    
%*******************************************************************************
% Target body name: Craft (NA)                      {source: DE431mx}
% Center body name: Earth (399)                     {source: DE431mx}
% Center-site name: BODY CENTER
%*******************************************************************************
ecend_c_leo2=0; %Eccentricity, e for circular orbits is zero
OM_c_leo2=0;%Longitude of Ascending Node, OM(degrees)    
w_c_leo2= 0; %Argument of Perifocus, w (degrees)
%-------------------------Orbital Elements Craft---------------------------

%First of we need to calculate the position of the craft in leo 1
[ri_c_leo1,~]=COE2RV(a_c_leo1,ecend_c_leo1,incl_c_leo1,OM_c_leo1,w_c_leo1,nu_c_leo1,mu_e);
%Since the craft vectors are geocentric, we need to transform them to
%Heleocentric. This can be done by simply adding the Heleocentric position
%vector of earth. So the earth's vector are

[ri_e,vi_e]=COE2RV(a_e,ecen_e,incl_e,OM_e,w_e,nu_e,mu_sun);
dv1=vi_e';
%Transforming craft vector by directly adding it to the earth's vector.
ri_e_bckp=ri_e;
ri_e=ri_e+ri_c_leo1; %Transform Earth position
%A crude assumtion above is that the heleocentric velicites of craft are
%negligeable compared to earth's.

[ri_d,vi_d]=COE2RV(a_d,ecen_d,incl_d,OM_d,w_d,nu_d,mu_sun);

%Solve Lamberts problem for the velocities requited to reach from Leo1 to
%EROS in dt0 time
[V1, V2, ~, ~] = lambert(ri_e', ri_d', dt0, 0, mu_sun);
C3_e2leo1=(norm(V2-V1))^2; %Characteristic energy
dv2=V1;dv3=V2;
%Now we need the orbital elements and ture anomalies for the orbit.
%Calculate initial transfer orbital elements
[ a,ecen,incl,Omega,w,nu0 ] = RV2COE(ri_e', V1, mu_sun ); 
%Calculate the true anomaly end
[ ~,~,~,~,~,nu01 ] = RV2COE(ri_d', V2, mu_sun );

if nu0<nu01
    [ris,~]=COE2RV(a,ecen,incl,Omega,w,nu01,mu_sun);
else
    [ris,~]=COE2RV(a,ecen,incl,Omega,w,nu0,mu_sun);
end
%Reset Earth position to orignal (i.e. remove the tranformation done before)
rie_e=ri_e_bckp;
%Let do the opposite
%First we need to calculate the leo2 earth insertion orbit.
[ri_c_leo2,~]=COE2RV(a_c_leo2,ecend_c_leo2,incl_c_leo2,OM_c_leo2,w_c_leo2,nu_c_leo2,mu_e);
%Now lets calculate the new positions after dt0
%Following is current position of eros. We only need to find the new state 
%vector after dt1
rd_ini=ris'; %Position vector where the statellite landed on EROS.
%Now we need to calculate the position of EROS (Thus satellite) after dt1
%amount of days.
[rd_ini,vd_ini] = rv_from_r0v0(rd_ini,vi_d', dt1*86400,mu_sun);
re_ini=ri_e; %This is the position when the satellite left LEO 1.
%Since the position of Earth has not been updated since dto,
%we need to calculate its position since dt0+dt1.
[re_ini,dv4] = rv_from_r0v0(re_ini,vi_e, (dt0+dt1)*86400,mu_sun);
dv4=dv4';
re_ini=re_ini+ri_c_leo2;
%Now we calculate the velocities required for transiting from EROS to Earth
[V1, V2, ~, ~] = lambert(rd_ini, re_ini', dt1, 0, mu_sun);

C3_e2leo2=(norm(V2-V1))^2; %Characteristic energy
C3_total=C3_e2leo1+C3_e2leo2;
dv5=V1;dv6=V2;

total=norm(dv6-dv5)+norm(dv3-dv2);
norm(dv2-dv1)+norm(dv3-dv2)+norm(dv4-dv3)+norm(dv5-dv4)+norm(dv6-dv5);
vector=(dv2-dv1)+(dv3-dv2)+(dv4-dv3)+(dv5-dv4)+(dv6-dv5);
end