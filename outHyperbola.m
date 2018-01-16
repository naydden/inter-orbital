clc;
clear all;

%% DATA
load('DATA.mat');
% geocentric velocity when going out of Earth Influence zone in the K system
Vinf = 3e5*[10,1,0];
% geocentric velocity in ECI
Vinf_eci = rotx(-earth.epsilon)*Vinf';
% normalized Vinf_eci
vinf_eci = Vinf_eci/norm(Vinf_eci);

%% parkingOrbit
h = 1000;
ro = earth.R + h;
Vo = sqrt(earth.mu/ro); % velocity at parking orbit;

%% deltaV
deltaV = sqrt(norm(Vinf)^2+2*Vo)-Vo;

%% hyperbolic path
a = mu/(Vu^2);
e = 1 + (Vinf/Vo)^2;
beta = acos(1/e);
b = a*sqrt(e^2-1);
% y = ((x+a*e)/(a^2)-1)*b^2;

%% elements
% r i v de la sonda on s'ha d'aplicar delta V. Punt I en figura 5.26.

%!! ARA R ES DEL PUNT C EN L_ORBITA. NO SE COM PASSAR A PUNT I (amb Beta)
r = ro*[ -vinf_eci(1), -vinf_eci(2), -vinf_eci(3)];
v = 1000*[-12,3,-2];
% element orbitals del punt I en sistema ECI
elem = rToElementsECI(r,v,earth);