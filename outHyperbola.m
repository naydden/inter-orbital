function [hyperbolaExit, deltaV] = outHyperbola (v_inf)
%  function that gets Vinf and gives hyperbolic trajectory
%  and the necessary deltaV

%% DATA
R_e = 6.3782e+03;  %km
mu_sun = 1.3271741784e20;
mu_e = 3.9820e+14; %SI
%% parkingOrbit
h = 800000; %height in m
ro = R_e*1000 + h;
Vo = sqrt(mu_e/ro); % velocity at parking orbit;
%% deltaV
Vinf = norm(v_inf);
deltaV = sqrt(Vinf^2+2*Vo^2)-Vo;
%% hyperbolic path
hyperbolaExit.a = mu_e/(Vinf^2);
hyperbolaExit.e = 1 + (Vinf/Vo)^2;
hyperbolaExit.beta = acosd(1/hyperbolaExit.e);
hyperbolaExit.b = hyperbolaExit.a*sqrt(hyperbolaExit.e^2-1);
end