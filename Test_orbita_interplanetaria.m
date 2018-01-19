% % Cas de la Terra a Mart
r1 = [0.4537 -0.9094 0];
r2 = [0.3148 1.5078 0.0239];
t1 = JulianDate(2020,7,19);
t2 = JulianDate(2021,1,25);
% 
% % Cas de Mart a Júpiter
% r1 = [1.3277 0.4901 0.0223];
% r2 = [5.0135 2.1380 0.0505];
% t1 = JulianDate(2026,6,5);
% t2 = JulianDate(2029,4,25);
% 
% % Cas 1 de Mart a Júpiter
% r1 = [1.0707 0.9868 0.0055];
% r2 = [5.2210 1.4357 0.1109];
% t1 = JulianDate(2037,10,25);
% t2 = JulianDate(2039,10,15);
% 
% % Cas 2 de la Terra a Mart
% r1 = [0.9848 0.1338 0];
% r2 = [0.6797 1.2298 0.0424];
% t1 = JulianDate(2033,3,13);
% t2 = JulianDate(2033,8,5);
% 
% % Cas 3 de la Terra a Mart
% r1 = [0.5264 0.8316 0.0001];
% r2 = [0.0108 1.4542 0.0309];
% t1 = JulianDate(2031,1,23);
% t2 = JulianDate(2031,8,1);
% 
% % Cas 4 de la Terra a Mart
% r1 = [0.4342 0.9188 0.0001];
% r2 = [0.6775 1.3571 0.0118];
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
% r1 = [1.5831 0.3913 0.0306];
% r2 = [0.9123 0.4340 0];
% t1 = JulianDate(2033,1,18);
% t2 = JulianDate(2033,8,28);
% 
% % Cas 7 de Mart a la Terra
% r1 = [1.4166 0.8722 0.0530];
% r2 = [0.2345 0.9893 0.0001];
% t1 = JulianDate(2030,11,20);
% t2 = JulianDate(2031,7,6);
% 
% % Cas 8 de la Terra a Mart (hiperbòlica)
% r1 = [0.4383 0.8843 0];
% r2 = [-0.2082 -1.4582 -0.0255];
% t1 = JulianDate(2021,11,26);
% t2 = JulianDate(2022,2,19);
% 
% % Cas 9 de la Terra a Mart (hiperbòlica)
% r1 = [-0.4079 0.8950 0];
% r2 = [0.6393 -1.2542 -0.0420];
% t1 = JulianDate(2022,1,15);
% t2 = JulianDate(2022,4,20);

deltat = t2-t1;

[a, e, theta1, w, i, Omega] = orbita_interplanetaria(r1,r2,deltat)
