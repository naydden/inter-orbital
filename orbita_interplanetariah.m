function [a, e, theta1, w, i, Omega] = orbita_interplanetariah(r1,r2,deltat)
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
lambda1 = checkTangent(lambda1,r1(1),r1(2));
lambda2 = atan(r2(2)/r2(1)); % [rad]
lambda2 = checkTangent(lambda2,r2(1),r2(2));
beta1 = asin(r1(3)/norm(r1)); % [rad]
beta2 = asin(r2(3)/norm(r2)); % [rad]

% Càlcul d'increment d'angles
deltalambda = wrapTo2Pi(lambda2-lambda1); % [rad]
deltatheta = acos(sin(beta1)*sin(beta2)+cos(beta1)*cos(beta2)*cos(deltalambda)); % [rad]

% Resolució de les equacions
syms e a theta1;
eqn1 = (norm(r2)-norm(r1))/(norm(r1)*cos(theta1)-norm(r2)*cos(theta1+deltatheta))-e == 0;
eqn2 = norm(r1)*(1+e*cos(theta1))/(e^2-1)-a == 0;
eqn3 = 365.25*a^(3/2)*(e*sqrt(e^2-1)*sin(theta1+deltatheta)/(1+e*cos(theta1+deltatheta))-log(((tan((theta1+deltatheta)/2)+sqrt((e+1)/(e-1)))/(tan((theta1+deltatheta)/2)-sqrt((e+1)/(e-1)))))-e*sqrt(e^2-1)*sin(theta1)/(1+e*cos(theta1))+log(((tan(theta1/2)+sqrt((e+1)/(e-1)))/(tan(theta1/2)-sqrt((e+1)/(e-1))))))/(2*pi)-deltat == 0;
S = solve(eqn1,eqn2,eqn3);

e = double(S.e); % excentricitat
a = double(S.a); % semieix major
theta1 = double(S.theta1); % posició inicial en l'òrbita [rad]

% Correcció per si surt e negatiu
if e<0
    e = -e;
    theta1 = theta1+pi;
end
theta1 = wrapTo2Pi(theta1);

% Càlcul de la inclinació per trigonometria esfèrica
A = asin(cos(beta2)*sin(deltalambda)/sin(deltatheta)); % [rad]
i = acos(sin(A)*cos(beta1)); % [rad]
if i>pi/2
    i = i-pi;
end
l = asin(tan(beta1)/tan(i)); % [rad]
if beta1<0 && l>0
    l = -l;
end

sigma = atan(tan(beta1)/cos(A)); % [rad]
Omega = lambda1-l; % [rad]
w = 2*pi-(theta1-sigma); % [rad]

if(i<0)
    i = abs(i);
    Omega = Omega+pi;
    w = w+pi;
end

theta1 = rad2deg(theta1); % [º]
i = rad2deg(i); % [º]
Omega = rad2deg(wrapTo2Pi(Omega)); % [º]
w = rad2deg(wrapTo2Pi(w)); % [º]

end