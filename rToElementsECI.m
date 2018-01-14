%function elem = rToElementsECI(r,v,mu)
clc;
clear all;
G = 6.67e-11;
Mt = 5.97e24;
%mu = G * Mt;
mu = 3.986005e14;
elem = struct('a',0,'e',0,'i',0,'Omega',0,'omega',0,'theta',0,'M',0);
r = 1000*[ 10000, 4000, -4000];
v = 1000*[-12,3,-2];
%% angular moment L and total energy E
L = cross(r,v);
E = 0.5*norm(v)^2-mu/norm(r);
%% eccentricity and conical parameter p and semi-major axis
e = sqrt(1+(2*E*(norm(L)^2)/(mu^2)));
p = (norm(L)^2)/mu;
if e<1
    a = p/(abs(1-e^2));
else
    a = p/(abs(e^2-1));
end
%% mean anomaly
esinTheta = norm(L)/(mu*norm(r))*dot(r,v);
ecosTheta = p/norm(r)-1;
if esinTheta > 0 && ecosTheta >0
    % primer quadrant
    theta = acos(ecosTheta/e);
elseif esinTheta > 0 && ecosTheta <0
    % segon quadrant
    theta = acos(ecosTheta/e);
elseif esinTheta < 0 && ecosTheta <0
    % tercer quadrant
    theta = pi + (pi - acos(ecosTheta/e));
else
    % quart quadrant
    theta = 2*pi - acos(ecosTheta/e);
end

% eccentric annomaly
numE = tan(theta/2);
denE = sqrt((1+e)/(1-e));
Ehalf = atan(numE/denE);
Ehalf = checkTangent(Ehalf,numE,denE);
% mean anomaly (Kepler)
M = Ehalf*2-e*sin(Ehalf*2);

%% Euler angles
W = L/norm(L);
i = atan((sqrt(W(1)^2+W(2)^2)/W(3)));
if W(3) < 0
    i = i + pi;
end
Omega = atan(W(1)/(-W(2)));
numO = W(1);
denO = -W(2);
Omega = checkTangent(Omega,numO,denO);
numU = r(3)/sin(i);
denU = r(1)*cos(Omega)+r(2)*sin(Omega);
u = atan((numU)/(denU));
u = checkTangent(u,numU,denU);
omega = u - theta;
%% final
elem.a = a/1000;
elem.e = e;
elem.i = i/pi*180;
elem.Omega = Omega/pi*180;
elem.omega = omega/pi*180;
elem.theta = theta/pi*180;
elem.M = M;
%end