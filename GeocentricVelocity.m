clc;
clear all;

% Sun
mu=1.3275e20;

% Time
t=(20*365.25+201)*24*3600;

% Earth orbit
a_e=149598023000;
e_e=0.0167086;
I_e=0.00005;
RAAN_e=348.73936;
AP_e=114.20783;
M0_e=358.617;
t0_e=0;

% Earth vectors
[r_e,v_e,theta_e] = OrbitalVectors (t,mu,a_e,e_e,I_e,RAAN_e,AP_e,M0_e,t0_e);

% Spacecraft orbit
a_s=199068390561.7830;
e_s=0.23629;
I_s=1.435;
RAAN_s=296.424;
AP_s=0.470;
M0_s=359.7033431;
t0_s=t;

% Earth vectors
[r_s,v_s,theta_s] = OrbitalVectors (t,mu,a_s,e_s,I_s,RAAN_s,AP_s,M0_s,t0_s);

% Geocentric velocity
Vinf=v_s-v_e;