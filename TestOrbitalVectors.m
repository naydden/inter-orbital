year=2033;
month=8;
day=5;

t=JulianDate(year,month,day); %Time

mu=1.3271741784e20; %Gravitational constant x Mass of central body

a=227939200000; %Major semiaxis
e=0.0934; %Excentricity
I=1.850; %Inclination
RAAN=49.558; %Right ascension of the ascending node
AP=286.502; %Argument of the perigee
M0=19.3564; %Mean anomaly at reference time
t0=JulianDate(2000,1,1); %Reference time

[r,v,theta] = OrbitalVectors (t,mu,a,e,I,RAAN,AP,M0,t0);
r=r/149597870700
norm(r)