function [r,v,theta] = OrbitalVectors (t,mu,a,e,I,RAAN,AP,M0,t0)
%Function that computes the position vector (r), the velocity vector (v)
%and the true anomaly (theta) for a given orbital parameters
%
%Inputs
%t: Time at which to compute the outputs, in JD [days]
%mu: Gravitational constant multipied by the mass of the central body (G*M)
%    [N*m^2/kg^2]
%a: Major semiaxis [m]
%e: Excentricity [dimensionless]
%I: Inclination [º]
%RAAN: Right ascension of the ascending node [º]
%AP: Argument of the perigee [º]
%M0: Mean anomaly at a reference time [º]
%t0: Reference time in JD [days]
%
%Outputs
%r: Vector for the posicion of the orbiting object in the same system of
%   reference as the orbital parameters [m]
%r: Vector for the velocity of the orbiting object in the same system of
%   reference as the orbital parameters [m/s]
%theta: True anomaly [º]

I=deg2rad(I); %Inclination in rad
RAAN=deg2rad(RAAN); %Right ascension of the ascending node in rad
AP=deg2rad(AP); %Argument of the perigee in rad
M0=deg2rad(M0); %Mean anomaly at reference time t0 in rad

T=sqrt(4*pi^2*a^3/mu); %Period [s]
n=2*pi/(T/3600/24); %Mean motion [rad/s]
M=M0+n*(t-t0); %Mean anomaly [rad]
M=wrapTo2Pi(M); %Mean anomaly between 0 and 2pi
if M>pi
    M=M-2*pi; %Correction for the hyperbole equations
end

E0=M+e*sin(M); %Initial excentric anomaly
if e<1 %Elyptic case
    p=a*(1-e^2); %Conic parameter
    for i=1:20 %Newton-Rapson
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
    for i=1:20 %Newton Rapson
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