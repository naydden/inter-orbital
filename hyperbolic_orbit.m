function [e,a,theta] = hyperbolic_orbit(r1,r2,deltat,deltatheta)
% Function that computes the excentricity, the semimajor axis and the
% true anomaly of an hyperbolic orbit using an iterative algorithm

% OUTPUTS
% e: eccentricity
% a: semi-major axis [AU]
% theta: true anomaly in t1 [rad]

% INPUTS
% r1: heliocentric position of the probe in t2 [AU]
% r2: heliocentric position of the probe in t2 [AU]
% deltat: t2-t1 [days]
% deltatheta: increment of the true anomaly between t1 and t2 [rad]

resta = 1000;
d = 1e-6; % error
theta = pi/2-acos(dot(r1,r2)/(norm(r1)*norm(r2))); % initial value
if theta>deg2rad(-50)
    theta = theta-deg2rad(10);
end
eant = 1000;

while(resta>d)
    e = (norm(r2)-norm(r1))/(norm(r1)*cos(theta)-norm(r2)*cos(theta+deltatheta));
    a = norm(r1)*(1+e*cos(theta))/(e^2-1);
    delta = 365.25*a^(3/2)*(e*sqrt(e^2-1)*sin(theta+deltatheta)/(1+e*cos(theta+deltatheta))-log(abs((tan((theta+deltatheta)/2)+sqrt((e+1)/(e-1)))/(tan((theta+deltatheta)/2)-sqrt((e+1)/(e-1)))))-e*sqrt(e^2-1)*sin(theta)/(1+e*cos(theta))+log(abs((tan(theta/2)+sqrt((e+1)/(e-1)))/(tan(theta/2)-sqrt((e+1)/(e-1))))))/(2*pi);
    
    resta = abs(e-eant);
    eant = e;
    
    if resta>d
        theta = theta+(deltat-delta)/10000;
    end
end
end
