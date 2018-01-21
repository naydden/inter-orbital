function [hyperbolaExit, parkingOrbit, deltaV] = outHyperbola (t,a_s,e_s,i_s,Omega_s,omega_s, theta1,v1,v_E,v_inf)
%  function that gets orbital elements of the elliptic trajectory and
%  obtains deltaV, the parking Orbit before exit and the hyperbolic
%  trajectory

% Outputs in ECI system:
%  hyperbolaExit (struct) -> a, e i b
%  parkingOrbit (struct) -> a,e,i,Omega,omega,theta,M
%  deltaV necessary is in Ecliptic coordinate system

%% DATA
load('DATA.mat');

% Earth orbit
a_e=149598023000;
e_e=0.0167086;
I_e=0.00005;
RAAN_e=348.73936;
AP_e=114.20783;
M0_e=358.617;
t0_e=0;

%% Velocities

% geocentric velocity in ECI
Vinf_eci = rotx(-earth.epsilon)*v_inf';
VTierra_eci = rotx(-earth.epsilon)*v_E';

% normalized Vinf_eci
vinf_eci = Vinf_eci/norm(Vinf_eci);

%% parkingOrbit
h = 1000;
ro = earth.R + h;
Vo = sqrt(earth.mu/ro); % velocity at parking orbit;

%% deltaV
deltaV = sqrt(norm(Vinf)^2+2*Vo)-Vo;

%% hyperbolic path
Vinf = norm(v_inf);
hyperbolaExit.a = mu/(Vu^2);
hyperbolaExit.e = 1 + (Vinf/Vo)^2;
beta = acos(1/e);
hyperbolaExit.b = a*sqrt(e^2-1);
% y = ((x+a*e)/(a^2)-1)*b^2;

%% parking orbit

% POINT C
% r is an oposed vector to vinf_eci with the radius of the orbit
r = ro*[ -vinf_eci(1), -vinf_eci(2), -vinf_eci(3)];
% v. 
n_plano = cross(VTierra_eci,Vinf_eci);
v = cross(r,n_plano);
% velocitat de la sonda per la direccio de v
v = Vo*v/norm(v);

elem = rToElementsECI(r,v,earth);

% In PQW system
r_pqw = [ro,0,0];
% rotation of r_pqw to beta angle. Point I
r_pqw_beta = rotz(pi+beta)*r_pqw';
% PQW2IJK
r_ijk_beta = rotz(-elem.Omega)*rotx(-elem.i)*rotz(-elem.omega)*r_pqw_beta;
% IJK2ECI
r_eci_beta = rotx(-earth.epsilon)*r_ijk_beta';
% Velocity at I
v_eci_beta = cross(r_eci_beta,n_plano);
v_eci_beta = Vo*v/norm(v_eci_beta);
r_eci_beta
v_eci_beta
parkingOrbit = rToElementsECI(r_eci_beta,v_eci_beta,earth);
end