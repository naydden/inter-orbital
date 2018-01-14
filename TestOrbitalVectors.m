t=0; %Time

mu=3.986005e14; %Gravitational constant x Mass of central body

a=25512556 ; %Major semiaxis
e=0.5; %Excentricity
I=30; %Inclination
RAAN=180; %Right ascension of the ascending node
AP=270; %Argument of the perigee
M0=0; %Mean anomaly at reference time
t0=0; %Reference time

[r,v,theta] = OrbitalVectors (t,mu,a,e,I,RAAN,AP,M0,t0);