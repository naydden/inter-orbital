% % Cas de la Terra a Mart
% r1 = [0.4537 -0.9094 0];
% r2 = [0.3148 1.5078 0.0239];
% t1 = JulianDate(2020,7,19);
% t2 = JulianDate(2021,1,25);
% 
% % Cas de Mart a J�piter
% r1 = [1.3277 0.4901 -0.0223];
% r2 = [-5.0135 -2.1380 0.0505];
% t1 = JulianDate(2026,6,5);
% t2 = JulianDate(2029,4,25);
% 
% Cas de la Terra a Mart (hiperb�lic)
r1 = [-0.9609 0.2466 0];
r2 = [0.7285 -1.1980 -0.0430];
t1 = JulianDate(2020,3,6);
t2 = JulianDate(2020,6,9);
% 
% % Cas 1 de Mart a J�piter
% r1 = [1.0707 0.9868 -0.0055];
% r2 = [-5.2210 1.4357 0.1109];
% t1 = JulianDate(2037,10,25);
% t2 = JulianDate(2039,10,15);
% 
% % Cas 2 de la Terra a Mart
% r1 = [-0.9848 0.1338 -0];
% r2 = [0.6797 -1.2298 -0.0424];
% t1 = JulianDate(2033,3,13);
% t2 = JulianDate(2033,8,5);
% 
% % Cas 3 de la Terra a Mart
% r1 = [-0.5264 0.8316 0.0001];
% r2 = [0.0108 -1.4542 -0.0309];
% t1 = JulianDate(2031,1,23);
% t2 = JulianDate(2031,8,1);
% 
% % Cas 4 de la Terra a Mart
% r1 = [0.4342 -0.9188 -0.0001];
% r2 = [-0.6775 -1.3571 -0.0118];
% t1 = JulianDate(2025,7,18);
% t2 = JulianDate(2025,10,21);
% 
% % Cas 5 de la Terra a Venus
% r1 = [-0.4255 -0.9194 0];
% r2 = [0.0356 0.7189 0.0079];
% t1 = JulianDate(2023,5,27);
% t2 = JulianDate(2023,11,1);
% 
% % Cas 6 de Mart a la Terra
% r1 = [-1.5831 -0.3913 0.0306];
% r2 = [0.9123 -0.4340 -0];
% t1 = JulianDate(2033,1,18);
% t2 = JulianDate(2033,8,28);
% 
% % Cas 7 de Mart a la Terra
% r1 = [-1.4166 0.8722 0.0530];
% r2 = [0.2345 -0.9893 -0.0001];
% t1 = JulianDate(2030,11,20);
% t2 = JulianDate(2031,7,6);
% 
% % Cas 8 de la Terra a Mart (hiperb�lica)
% r1 = [0.4383 0.8843 0];
% r2 = [-0.2082 -1.4582 -0.0255];
% t1 = JulianDate(2021,11,26);
% t2 = JulianDate(2022,2,19);
% 
% % Cas 9 de la Terra a Mart (hiperb�lica)
% r1 = [-0.4079 0.8950 0];
% r2 = [0.6393 -1.2542 -0.0420];
% t1 = JulianDate(2022,1,15);
% t2 = JulianDate(2022,4,20);

deltat = t2-t1;

[a, e, theta1, w, i, Omega] = orbita_interplanetaria(r1,r2,deltat)

lambda1 = atan(r1(2)/r1(1)); %  [rad]
lambda1 = checkTangent(lambda1,r1(2),r1(1));
lambda2 = atan(r2(2)/r2(1)); % [rad]
lambda2 = checkTangent(lambda2,r2(2),r2(1));
beta1 = asin(r1(3)/norm(r1)); % [rad]
beta2 = asin(r2(3)/norm(r2)); % [rad]
deltalambda = wrapTo2Pi(lambda2-lambda1); % [rad]
deltatheta = acos(sin(beta1)*sin(beta2)+cos(beta1)*cos(beta2)*cos(deltalambda)); % [rad]

% Mart
a_M = 227939200000/149597870700;
e_M = 0.0934;
t_M=0:pi/40:2*pi;
R_M=a_M*(1-e_M^2)./(1+e_M*cos(t_M));
polar(t_M,R_M);
hold on

% Terra
a_E = 149598023000/149597870700;
e_E = 0.0167086;
t_E=0:pi/40:2*pi;
R_E=a_E*(1-e_E^2)./(1+e_E*cos(t_E));
polar(t_E,R_E);
hold on

% Interplanet�ria
t=deg2rad(theta1):pi/40:deg2rad(theta1)+deltatheta;
if e<1
    R=a*(1-e^2)./(1+e*cos(t));
else
    R=a*(e^2-1)./(1+e*cos(t));
end
polar(t,R);

legend('Mars','Earth','Interplanetary')