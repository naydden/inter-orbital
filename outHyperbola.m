%% DATA
Rt = 6371;
G = 6.67e-11;
Mt = 5.97e24;
mu = G * Mt;
Vinf = 3e5*[10,1,0]; % geocentric velocity when going out of Earth Influence zone

%% parkingOrbit
h = 1000;
r = Rt + h;
Vo = sqrt(mu/r); % velocity at parking orbit;

a = r;
e = 0;
i = 
Omega = 
omega = 

%% deltaV
deltaV = sqrt(norm(Vinf)^2+2*Vo)-Vo;

%% hyperbolic path
a = mu/(Vu^2);
e = 1 + (Vinf/Vo)^2;
beta = acos(1/e);
b = a*sqrt(e^2-1);
% y = ((x+a*e)/(a^2)-1)*b^2;