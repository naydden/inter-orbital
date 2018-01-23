function [r,v,theta] = OrbitalVectors (t,mu,a,e,I,RAAN,AP,M0,t0)
% Function that computes the position vector (r), the velocity vector
% (v) and the true anomaly (theta) for a given orbital parameters

% OUTPUTS
% r: position vector of the orbiting object in the same system of
% reference as the orbital parameters [m]
% v: velocity vectorof the orbiting object in the same system of
% reference as the orbital parameters [m/s]
% theta: true anomaly [deg]

% INPUTS
% t: time at which to compute the outputs in JD [days]
% mu: gravitational constant multiplied by the mass of the central body
% (G*M) [N*m^2/kg^2]
% a: semi-major axis [m]
% e: eccentricity
% I: inclination [deg]
% RAAN: right ascension of the ascending node [deg]
% AP: argument of the perigee [deg]
% M0: mean anomaly at a reference time [deg]
% t0: reference time in JD [days]

I=deg2rad(I); %Inclination [rad]
RAAN=deg2rad(RAAN); %Right ascension of the ascending node [rad]
AP=deg2rad(AP); %Argument of the perigee [rad]
M0=deg2rad(M0); %Mean anomaly at reference time t0 [rad]
t=t*24*3600; %Time [s]
t0=t0*24*3600; %Reference time [s]

T=sqrt(4*pi^2*a^3/mu); %Period [s]
n=2*pi/T; %Mean motion [rad/s]
M=M0+n*(t-t0); %Mean anomaly [rad]
M=wrapTo2Pi(M); %Mean anomaly between 0 and 2pi
if M>pi
    M=M-2*pi; %Correction for the hyperbolic equations
end

E0=M+e*sin(M); %Initial eccentric anomaly
error=1e-8;
if e<1 %Elliptic case
    p=a*(1-e^2); %Conic parameter
    E=1;
    while abs(E-E0)>error %Newton-Rapson
        E=E0+(M-E0+e*sin(E0))/(1-e*sin(E0));
        E0=E;
    end
    theta=2*atan(sqrt((1+e)/(1-e))*tan(E/2)); %True anomaly
    r_mod=a*(1-e*cos(E)); %Modulus of the position vector
else %Hyperbolic case
    p=a*(e^2-1); %Conic parameter
    if e<1.6     
        if M<=pi
            F0=M+e; %Initial hyperbolic anomaly
        else
            F0=M-e; %Initial hyperbolic anomaly
        end
    else
        if e<3.6 && M>pi
            F0=M-e; %Initial hyperbolic anomaly
        else
            F0=M/(e-1); %Initial hyperbolic anomaly
        end
    end
    F=1;
    while abs(F-F0)>error %Newton Rapson
        F=F0+(M-e*sinh(F0)+F0)/(e*cosh(F0)-1);
        F0=F;
    end
    theta=2*atan(sqrt((e+1)/(e-1))*tanh(F/2)); %True anomaly
    r_mod=a*(e*cosh(F)-1); %Modulus of the position vector
end

%Rotation coefficients
Px=cos(RAAN)*cos(AP)-sin(RAAN)*cos(I)*sin(AP);
Py=sin(RAAN)*cos(AP)+cos(RAAN)*cos(I)*sin(AP);
Pz=sin(I)*sin(AP);
Qx=-cos(RAAN)*sin(AP)-sin(RAAN)*cos(I)*cos(AP);
Qy=-sin(RAAN)*sin(AP)+cos(RAAN)*cos(I)*cos(AP);
Qz=sin(I)*cos(AP);
Wx=sin(RAAN)*sin(I);
Wy=-cos(RAAN)*sin(I);
Wz=cos(I);
P=[Px ; Py ; Pz]; %Rotation vector for x_orb
Q=[Qx ; Qy ; Qz]; %Rotation vector for y_orb

r=r_mod*(cos(theta)*P+sin(theta)*Q); %Distance vector
v=sqrt(mu/p)*(-sin(theta)*P+(e+cos(theta))*Q); %velocity vector
end