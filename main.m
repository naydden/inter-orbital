
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

r1 = [0.4537 -0.9094 0];
r2 = [0.3148 1.5078 0.0239];
t1 = JulianDate(2020,7,19);
t2 = JulianDate(2021,1,25);

deltat = t2-t1;

%% Part 1: Heliocentric elliptic trajectory
[a, e, theta1, w, i, Omega] = orbita_interplanetaria(r1,r2,deltat);

%% Part 2: Exit. Parking orbit, deltaV and hyperbolic trajectory

[hyperbolaExit, parkingOrbit, deltaV] = outHyperbola (t1,a,e,i,Omega,w,theta1);

%% Part 3: Arrival. Hyperbolic trajectory, deltaV and parking orbit

%% Part 4: Results Presentation