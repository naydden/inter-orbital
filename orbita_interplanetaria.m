function [a, e, theta1, w, i, Omega] = orbita_interplanetaria(r1,r2,deltat)
% a: semieix major [AU]
% e: excentricitat
% theta1: angle entre el punt Aries i la sonda (o planeta de sortida) [º]
% w: argument del periheli [º]
% i: inclinació [º]
% Omega: longitud eclíptica del node ascendent [º]

% r1: posició heliocèntrica de la sonda en t1 [AU]
% r2: posició heliocèntrica de la sonda en t2 [AU]
% deltat: t2-t1 [dies]

% Càlcul d'angles
lambda1 = atan(r1(2)/r1(1)); %  [rad]
% lambda1 = checkTangent(lambda1,r1(2),r1(1));
lambda2 = atan(r2(2)/r2(1)); % [rad]
% lambda2 = checkTangent(lambda2,r2(2),r2(1));
beta1 = asin(r1(3)/norm(r1)); % [rad]
beta2 = asin(r2(3)/norm(r2)); % [rad]

% Càlcul d'increment d'angles
deltalambda = lambda2-lambda1; % [rad]
deltatheta = acos(sin(beta1)*sin(beta2)+cos(beta1)*cos(beta2)*cos(deltalambda)); % [rad]

% Resolució de les equacions
syms e a theta1;
eqn1 = (norm(r2)-norm(r1))/(norm(r1)*cos(theta1)-norm(r2)*cos(theta1+deltatheta))-e == 0;
eqn2 = norm(r1)*(1+e*cos(theta1))/(1-e^2)-a == 0;
eqn3 = 365.25*a^(3/2)*(2*atan(sqrt((1-e)/(1+e))*tan((theta1+deltatheta)/2))-e*sqrt(1-e^2)*sin(theta1+deltatheta)/(1+e*cos(theta1+deltatheta))-2*atan(sqrt((1-e)/(1+e))*tan(theta1/2))+e*sqrt(1-e^2)*sin(theta1)/(1+e*cos(theta1)))/(2*pi)-deltat == 0;
S = solve(eqn1,eqn2,eqn3);

e = S.e; % excentricitat
a = S.a; % semieix major
theta1 = wrapTo360(rad2deg(S.theta1)); % posició inicial en l'òrbita [º]

% Càlcul de la inclinació per trigonometria esfèrica
A = asin(cos(beta2)*sin(deltalambda)/sin(deltatheta)); % [rad]
i = acos(sin(A)*cos(beta1)); % [rad]
l = asin(tan(beta1)/tan(i)); % [rad]

%  we need to check the members of the tangent here in order to decide
%  sigma
sigma = atan(tan(beta1)/cos(A)); % [rad]

Omega = lambda1-l; % [rad]
w = 2*pi-(deg2rad(theta1)-sigma); % [rad]

if(i<0)
    i = abs(i);
    Omega = Omega+pi;
    w = w+pi;
end

i = rad2deg(i); % [º]
Omega = wrapTo360(rad2deg(Omega)); % [º]
w = wrapTo360(rad2deg(w)); % [º]

end