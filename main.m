
%% Patched Conics. Interplanetary Trajectories.
% Authors:
%   - Silvia Gonzalez
%   - Laura Pla
%   - Josep Maria Serra Moncunill
%   - Boyan Naydenov
% Subject: Aerodynamics, Flight and Orbital Mechanics.
% Date: January 20th, 2018

%%
clear
clc

%% Data
load('DATA.mat');
mu = 1.3271741784e20;
dST = 149597870700;
%Earth data
a_E = 149598023000;
e_E = 0.0167086;
I_E = 0.00005;
RAAN_E = -11.26064;
AP_E = 114.20783;
M0_E = 358.617;
a_O = 149598023000;

%Mars data
a_M = 227939200000;
e_M = 0.0934;
I_M = 1.850;
RAAN_M = 49.558;
AP_M = 286.502;
M0_M = 19.3564;

% Time data
t0 = JulianDate(2000,1,1);
t1 = JulianDate(2020,7,19);
t2 = JulianDate(2021,1,25);
deltat = t2-t1;
T=sqrt(4*pi^2*a_O^3/mu)/(3600*24); %Period [days]

% Distances, velocities and true anomalies
[r1_E,v1_E,theta1_E] = OrbitalVectors (t1,mu,a_E,e_E,I_E,RAAN_E,AP_E,M0_E,t0);
r1_E = r1_E/dST;
[r2_M,v2_M,theta2_M] = OrbitalVectors (t2,mu,a_M,e_M,I_M,RAAN_M,AP_M,M0_M,t0);
r2_M = r2_M/dST;

%% Part 1: Heliocentric elliptic trajectory
[a_S, e_S, theta1, AP_S, I_S, RAAN_S] = orbita_interplanetaria(r1_E,r2_M,deltat);
E_S=acosd((e_S+cosd(theta1))/(1+e_S*cosd(theta1)));
M_S=E_S-e_S*sin(E_S);
a_S = a_S*dST;

%% Part 2: Exit. Geocentric Parking orbit and hyperbolic trajectory and deltaV
[r1,v1,theta1] = OrbitalVectors (t1,mu,a_S,e_S,I_S,RAAN_S,AP_S,M_S,t1);
v_inf1=v1-v1_E;
[hyperbolaExit, deltaV] = outHyperbola (t1,a_S,e_S,I_S,RAAN_S,AP_S,theta1,v1,v1_E,v_inf1);

%% Part 3: Arrival. Geocentric hyperbolic trajectory and Parking orbit and deltaV
[r2,v2,theta2] = OrbitalVectors (t2,mu,a_S,e_S,I_S,RAAN_S,AP_S,M_S,t1);
v_inf2=v2-v2_M;

%% Part 4: Results Presentation

N=100;
A=2;
t=linspace(t1-A*deltat,t2+A*deltat,(2*A+1)*N);
r_E=zeros(1,(2*A+1)*N,3);
r_M=zeros(1,(2*A+1)*N,3);
r_S=zeros(1,N,3);
for i=1:(2*A+1)*N
    [r_E(1,i,:),v1_E,theta1_E] = OrbitalVectors (t(i),mu,a_E,e_E,I_E,RAAN_E,AP_E,M0_E,t0);
    [r_M(1,i,:),v1_M,theta1_M] = OrbitalVectors (t(i),mu,a_M,e_M,I_M,RAAN_M,AP_M,M0_M,t0);
    if t(i)>=t1 && t(i)<=t2
        [r_S(1,i-A*N,:),v1_S,theta1_S] = OrbitalVectors (t(i),mu,a_S,e_S,I_S,RAAN_S,AP_S,M_S,t1);
    end
end
plot3(r_E(1,:,1),r_E(1,:,2),r_E(1,:,3),':b');
hold on
plot3(r_E(1,A*N:(A+1)*N,1),r_E(1,A*N:(A+1)*N,2),r_E(1,A*N:(A+1)*N,3),'-b','LineWidth',2);
hold on
plot3(r_M(1,:,1),r_M(1,:,2),r_M(1,:,3),':r');
hold on
plot3(r_M(1,A*N:(A+1)*N,1),r_M(1,A*N:(A+1)*N,2),r_M(1,A*N:(A+1)*N,3),'-r','LineWidth',2);
hold on
plot3(r_S(1,:,1),r_S(1,:,2),r_S(1,:,3),'-g','LineWidth',2);
axis equal;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.ZAxisLocation = 'origin';