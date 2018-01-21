
%% Patched Conics. Interplanetary Trajectories.
% Authors:
%   - Silvia Gonzalez
%   - Laura Pla
%   - Josep Maria Serra Moncunill
%   - Boyan Naydenov
% Subject: Aerodynamics, Flight and Orbital Mechanics.
% Date: January 20th, 2018

%%
clear
clc

%% Data
load('DATA.mat');
mu = 1.3271741784e20;

%Earth data
a_E = 149598023000;
e_E = 0.0167086;
I_E = 0.00005;
RAAN_E = -11.26064;
AP_E = 114.20783;
M0_E = 358.617;
a_O = 149598023000;

%Mars data
a_M = 227939200000;
e_M = 0.0934;
I_M = 1.850;
RAAN_M = 49.558;
AP_M = 286.502;
M0_M = 19.3564;

% Time data
t0 = JulianDate(2000,1,1);
t1 = JulianDate(2020,7,19);
t2 = JulianDate(2021,1,25);
deltat = t2-t1;
T=sqrt(4*pi^2*a_O^3/mu)/(3600*24); %Period [days]

% Distances, velocities and true anomalies
[r1_E,v1_E,theta1_E] = OrbitalVectors (t1,mu,a_E,e_E,I_E,RAAN_E,AP_E,M0_E,t0);
r1_E = r1_E/149597870700
[r2_M,v2_M,theta2_M] = OrbitalVectors (t2,mu,a_M,e_M,I_M,RAAN_M,AP_M,M0_M,t0);
r2_M = r2_M/149597870700

%% Part 1: Heliocentric elliptic trajectory
[a, e, theta1, w, i, Omega] = orbita_interplanetaria(r1_E,r2_M,deltat,T)
E=acosd((e+cosd(theta1))/(1+e*cosd(theta1)));
M=E-e*sin(E);

%% Part 2: Exit. Geocentric Parking orbit and hyperbolic trajectory and deltaV
[r1,v1,theta1] = OrbitalVectors (t1,mu,a,e,I,RAAN,AP,M,t1)
v_inf1=v1-v1_E;
[hyperbolaExit, parkingOrbit, deltaV] = outHyperbola (t1,a,e,i,Omega,w,theta1,v1,v1_E,v_inf1)

%% Part 3: Arrival. Geocentric hyperbolic trajectory and Parking orbit and deltaV
[r2,v2,theta2] = OrbitalVectors (t2,mu,a,e,I,RAAN,AP,M,t1);
v_inf2=v2-v2_M;

%% Part 4: Results Presentation