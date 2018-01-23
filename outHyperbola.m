function [hyperbolaExit, deltaV] = outHyperbola (v_inf)
%  Function that gets Vinf and gives hyperbolic trajectory
%  and the necessary deltaV

% OUTPUTS
% hyperbolaExit: orbital parameters of the hyperbola
% deltaV: Increment of velocity required to go from the parking orbit to
% the hyperbolic orbit [m/s]

% INPUT
% v_inf: planetocentric  hyperbolic excess velocity [m/s]

%% DATA
R_e = 6.3782e+03;  % [km]
mu_sun = 1.3271741784e20;
mu_e = 3.9820e+14; %SI

%% parkingOrbit
h = 800000; %height [m]
ro = R_e*1000 + h;
Vo = sqrt(mu_e/ro); % velocity in the parking orbit

%% deltaV
Vinf = norm(v_inf);
deltaV = sqrt(Vinf^2+2*Vo^2)-Vo;

%% hyperbolic path
hyperbolaExit.a = mu_e/(Vinf^2); % semi-major axis
hyperbolaExit.e = 1 + (Vinf/Vo)^2; % eccentricity
hyperbolaExit.beta = acosd(1/hyperbolaExit.e); % hyperbolic angle
hyperbolaExit.b = hyperbolaExit.a*sqrt(hyperbolaExit.e^2-1); % Impact parameter (semi-minor axis)
end