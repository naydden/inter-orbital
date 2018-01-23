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
mu = 1.3271741784e20;
dST = 149597870700;

%Venus data
a_V = 0.723330*dST;
e_V = 0.006772;
I_V = 3.3947;
RAAN_V = 76.6799;
AP_V = 54.8838;
M0_V = 50.4161;

%Earth data
a_E = 149598023000;
e_E = 0.0167086;
I_E = 0.00005;
RAAN_E = -11.26064;
AP_E = 114.20783;
M0_E = 358.617;

%Mars data
a_M = 227939200000;
e_M = 0.0934;
I_M = 1.850;
RAAN_M = 49.558;
AP_M = 286.502;
M0_M = 19.3564;

%Jupiter data
a_J = 7.7830e+11;
e_J = 0.048498;
I_J = 1.3033;
RAAN_J = 100.46;
AP_J = 273.866;
M0_J = 20.021;

% Time data
t0 = JulianDate(2000,1,1);
t1 = JulianDate(2020,7,19);
t2 = JulianDate(2021,1,25);
deltat = t2-t1;

% Distances, velocities and true anomalies (O:origin, D:destination)
[r1_O,v1_O,theta1_O] = OrbitalVectors (t1,mu,a_E,e_E,I_E,RAAN_E,AP_E,M0_E,t0);
r1_O = r1_O/dST;
[r2_D,v2_D,theta2_D] = OrbitalVectors (t2,mu,a_M,e_M,I_M,RAAN_M,AP_M,M0_M,t0);
r2_D = r2_D/dST;

%% Part 1: Heliocentric elliptic trajectory
[a_S, e_S, theta1, AP_S, I_S, RAAN_S] = orbita_interplanetaria(r1_O,r2_D,deltat);
E_S=acosd((e_S+cosd(theta1))/(1+e_S*cosd(theta1)));
M_S=E_S-e_S*sin(E_S);
a_S = a_S*dST;

%% Part 2: Exit. Geocentric Parking orbit and hyperbolic trajectory and deltaV
[r1,v1,theta1] = OrbitalVectors (t1,mu,a_S,e_S,I_S,RAAN_S,AP_S,M_S,t1);
v_inf1=v1-v1_O;
[hyperbolaExit, deltaV] = outHyperbola (v_inf1);

%% Part 3: Arrival. Geocentric hyperbolic trajectory and Parking orbit and deltaV
[r2,v2,theta2] = OrbitalVectors (t2,mu,a_S,e_S,I_S,RAAN_S,AP_S,M_S,t1);
v_inf2=v2-v2_D;

%% Part 4: Results Presentation
N=100;
A=2;
t=linspace(t1-A*deltat,t2+A*deltat,(2*A+1)*N);
r_O=zeros(1,(2*A+1)*N,3);
r_D=zeros(1,(2*A+1)*N,3);
r_S=zeros(1,N,3);
for i=1:(2*A+1)*N
    [r_O(1,i,:),v_O,theta_O] = OrbitalVectors (t(i),mu,a_E,e_E,I_E,RAAN_E,AP_E,M0_E,t0);
    [r_D(1,i,:),v_D,theta_D] = OrbitalVectors (t(i),mu,a_M,e_M,I_M,RAAN_M,AP_M,M0_M,t0);
    if t(i)>=t1 && t(i)<=t2
        [r_S(1,i-A*N,:),v1_S,theta1_S] = OrbitalVectors (t(i),mu,a_S,e_S,I_S,RAAN_S,AP_S,M_S,t1);
    end
end
plot3(r_O(1,:,1),r_O(1,:,2),r_O(1,:,3),':b');
hold on
plot3(r_O(1,A*N:(A+1)*N,1),r_O(1,A*N:(A+1)*N,2),r_O(1,A*N:(A+1)*N,3),'-b','LineWidth',2);
hold on
plot3(r_D(1,:,1),r_D(1,:,2),r_D(1,:,3),':r');
hold on
plot3(r_D(1,A*N:(A+1)*N,1),r_D(1,A*N:(A+1)*N,2),r_D(1,A*N:(A+1)*N,3),'-r','LineWidth',2);
hold on
plot3(r_S(1,:,1),r_S(1,:,2),r_S(1,:,3),'-g','LineWidth',2);
axis equal;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.ZAxisLocation = 'origin';