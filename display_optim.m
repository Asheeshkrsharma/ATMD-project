function display_optim(chromosome)
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

%**********************Mission Time frame (Since epoch)********************

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
OM_c_leo1=0;%Longitude of Ascending Node, OM(degrees)    
w_c_leo1= 0; %Argument of Perifocus, w (degrees)
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
orbit=0:0.07:2.2*pi; %True anomaly range is 0 to 360 degrees
orb_leo1=zeros(length(orbit),3);
%Lets also obtain the orbit.
for i=1:length(orbit)
    [orb_leo1(i,:),~]=COE2RV(a_c_leo1,ecend_c_leo1,incl_c_leo1,OM_c_leo1,w_c_leo1,orbit(i),mu_e);
end
%Since the craft vectors are geocentric, we need to transform them to
%Heleocentric. This can be done by simply adding the Heleocentric position
%vector of earth. So the earth's vector are

[ri_e,vi_e]=COE2RV(a_e,ecen_e,incl_e,OM_e,w_e,nu_e,mu_sun);

%Transforming craft vector by directly adding it to the earth's vector.
orb_leo1=bsxfun(@plus,orb_leo1,ri_e'); %Transform orbit
ri_e_bckp=ri_e;
ri_e=ri_e+ri_c_leo1; %Transform Earth position
%A crude assumtion above is that the heleocentric velicites of craft are
%negligeable compared to earth's.

[ri_d,vi_d]=COE2RV(a_d,ecen_d,incl_d,OM_d,w_d,nu_d,mu_sun);

%Solve Lamberts problem for the velocities requited to reach from Leo1 to
%EROS in dt0 time
[V1, V2, ~, ~] = lambert(ri_e', ri_d', dt0, 0, mu_sun);
C3=(norm(V1'-vi_e))^2;
vinf=norm(V2'-vi_d);

%Now we need the orbital elements and ture anomalies for the orbit.
%Calculate initial transfer orbital elements
[ a,ecen,incl,Omega,w,nu0 ] = RV2COE(ri_e', V1, mu_sun ); 
%Calculate the true anomaly end
[ ~,~,~,~,~,nu01 ] = RV2COE(ri_d', V2, mu_sun );

if nu0<nu01
    orbit=nu0:0.001:nu01;  %This is the range of true anomalies for the established orbit
else
    orbit=nu01:0.001:nu0;  %This is the range of true anomalies for the established orbit
end
ris=zeros(length(orbit),3); %cartesian orbit of satelite initialisation
ri_2e=zeros(length(orbit),3); %cartesian orbit of earth initialisation
ri_d=zeros(length(orbit),3); %cartesian orbit of EROS intialisation
scale=length(orbit)/360; %scale factor for earth and eros orbits
%Calculate the orbit.
for i=1:length(orbit)
   [ris(i,:),~]=COE2RV(a,ecen,incl,Omega,w,orbit(i),mu_sun);
   [ri_2e(i,:),~]=COE2RV(a_e,ecen_e,incl_e,OM_e,w_e,deg2rad(i+scale),mu_sun);
   [ri_d(i,:),~]=COE2RV(a_d,ecen_d,incl_d,OM_d,w_d,deg2rad(i+scale),mu_sun);   
end

%Plot transit orbit from LEO1 to EROS
hold on;
axis equal;
set(gca,'color','black')
set(gcf,'color','black')
set(gca,'XColor','white')
set(gca,'YColor','white')
set(gca,'ZColor','white')
[X Y Z] = sphere(50);          % Reference Sphere
Rm = 20*695700*ones(1,3);  % Mean Radius
HSUN = surf(Rm(1)*X, Rm(1)*Y, Rm(1)*Z);
        topoSUN = imread('euvisdoCarringtonMap.jpg');
        % Set it on SUN
        set(HSUN,'facecolor','texture',...
                'cdata',im2double(topoSUN),...
                'edgecolor','none');
        
plot3(ris(1,1),ris(1,2),ris(1,3),'w*','LineWidth',2,'MarkerSize',2)
plot3(ris(:,1),ris(:,2),ris(:,3),'g')
plot3(ris(1,1),ris(1,2),ris(1,3),'w*','LineWidth',2,'MarkerSize',2)
plot3(ris(end,1),ris(end,2),ris(end,3),'w+','LineWidth',2,'MarkerSize',7)
plot3(ri_2e(:,1),ri_2e(:,2),ri_2e(:,3),'b')
plot3(ri_d(:,1),ri_d(:,2),ri_d(:,3),'w')
plot3(orb_leo1(:,1),orb_leo1(:,2),orb_leo1(:,3),'w')
Re=ones(1,3)*Re;
HEARTH = surf(ri_e_bckp(1)+Re(1)*X,...
                      ri_e_bckp(2)+Re(2)*Y,...
                      ri_e_bckp(3)+Re(3)*Z);
topoEarth = imread('blackmarble_2016_3km.jpg');
set(HEARTH,'facecolor','texture',...
                'cdata',im2double(topoEarth),...
                'edgecolor','none');

plot3(ri_e_bckp(1),ri_e_bckp(2),ri_e_bckp(3),'b*','MarkerSize',15)

%Reset Earth position to orignal (i.e. remove the tranformation done before)
rie_e=ri_e_bckp;
%Let do the opposite

%First we need to calculate the leo2 earth insertion orbit.

[ri_c_leo2,~]=COE2RV(a_c_leo2,ecend_c_leo2,incl_c_leo2,OM_c_leo2,w_c_leo2,nu_c_leo2,mu_e);

%Lets also obtain the orbit.
orbit=0:0.07:2.2*pi; %True anomaly range is 0 to 360 degrees
orb_leo2=zeros(length(orbit),3);
for i=1:length(orbit)
    [orb_leo2(i,:),~]=COE2RV(a_c_leo2,ecend_c_leo2,incl_c_leo2,OM_c_leo2,w_c_leo2,orbit(i),mu_e);
end

%Now lets calculate the new positions after dt0
%Following is current position of eros. We only need to find the new state 
%vector after dt1
rd_ini=ris(end,:); %Position vector where the statellite landed on EROS.

%Now we need to calculate the position of EROS (Thus satellite) after dt1
%amount of days.
[rd_ini,~] = rv_from_r0v0(rd_ini,vi_d', dt1*86400,mu_sun);

re_ini=ri_e; %This is the position when the satellite left LEO 1.
%Since the position of Earth has not been updated since dto,
%we need to calculate its position since dt0+dt1.
[re_ini,~] = rv_from_r0v0(re_ini,vi_e, (dt0+dt1)*86400,mu_sun);

orb_leo2=bsxfun(@plus,orb_leo2,re_ini'); %Transform orbit
re_ini=re_ini+ri_c_leo2;

%Now we calculate the velocities required for transiting from EROS to Earth
[V1, V2, ~, ~] = lambert(rd_ini, re_ini', dt1, 0, mu_sun);

%Now we get the orbital elements for the Transit
[a,ecen,incl,Omega,w,nu0 ] = RV2COE(rd_ini, V1, mu_sun );
[ ~,~,~,~,~,nu01 ] = RV2COE(re_ini', V2, mu_sun );


%Plot
if nu0<nu01
    orbit=nu0:0.001:nu01;  %This is the range of true anomalies for the established orbit
else
    orbit=nu01:0.001:nu0;  %This is the range of true anomalies for the established orbit
end
ris=zeros(length(orbit),3); %cartesian orbit of satelite initialisation

for i=1:length(orbit)
   [ris(i,:),~]=COE2RV(a,ecen,incl,Omega,w,orbit(i),mu_sun);
end

plot3(rd_ini(1),rd_ini(2),rd_ini(3),'g+','LineWidth',2,'MarkerSize',12)
% plot3(re_ini(1),re_ini(2),re_ini(3),'w+','LineWidth',2,'MarkerSize',2)
plot3(ris(:,1),ris(:,2),ris(:,3),'r')
plot3(ris(1,1),ris(1,2),ris(1,3),'b+','LineWidth',2,'MarkerSize',12)
plot3(orb_leo2(:,1),orb_leo2(:,2),orb_leo2(:,3),'y')

