year=2023;
month=5;
day=27;

t=JulianDate(year,month,day); %Time

mu=1.3271741784e20; %Gravitational constant x Mass of central body

a=149598023000; %Major semiaxis
e=0.0167086; %Excentricity
I=0.00005; %Inclination
RAAN=-11.26064; %Right ascension of the ascending node
AP=114.20783; %Argument of the perigee
M0=358.617; %Mean anomaly at reference time
t0=JulianDate(2000,1,1); %Reference time

[r,v,theta] = OrbitalVectors (t,mu,a,e,I,RAAN,AP,M0,t0);
r=r/149597870700
norm(r)