function [a, e, theta1, w, i, Omega] = orbita_interplanetaria(r1,r2,deltat)
% a: semieix major [AU]
% e: excentricitat
% theta1: angle entre el punt Aries i la sonda (o planeta de sortida) [�]
% w: argument del periheli [�]
% i: inclinaci� [�]
% Omega: longitud ecl�ptica del node ascendent [�]

% r1: posici� helioc�ntrica de la sonda en t1 [AU]
% r2: posici� helioc�ntrica de la sonda en t2 [AU]
% deltat: t2-t1 [dies]

% C�lcul d'angles
lambda1 = atan(r1(2)/r1(1)); %  [rad]
lambda1 = checkTangent(lambda1,r1(2),r1(1));
lambda2 = atan(r2(2)/r2(1)); % [rad]
lambda2 = checkTangent(lambda2,r2(2),r2(1));
beta1 = asin(r1(3)/norm(r1)); % [rad]
beta2 = asin(r2(3)/norm(r2)); % [rad]

% C�lcul d'increment d'angles
deltalambda = wrapTo2Pi(lambda2-lambda1); % [rad]
deltatheta = acos(sin(beta1)*sin(beta2)+cos(beta1)*cos(beta2)*cos(deltalambda)); % [rad]

% Resoluci� de les equacions CAS EL�L�PTIC
syms e a theta1;
eqn1 = (norm(r2)-norm(r1))/(norm(r1)*cos(theta1)-norm(r2)*cos(theta1+deltatheta))-e == 0;
eqn2 = norm(r1)*(1+e*cos(theta1))/(1-e^2)-a == 0;
eqn3 = 365.25*a^(3/2)*(2*atan(sqrt((1-e)/(1+e))*tan((theta1+deltatheta)/2))-e*sqrt(1-e^2)*sin(theta1+deltatheta)/(1+e*cos(theta1+deltatheta))-2*atan(sqrt((1-e)/(1+e))*tan(theta1/2))+e*sqrt(1-e^2)*sin(theta1)/(1+e*cos(theta1)))/(2*pi)-deltat == 0;
S = solve(eqn1,eqn2,eqn3);

e = double(S.e); % excentricitat
a = double(S.a); % semieix major
theta1 = double(S.theta1); % posici� inicial en l'�rbita [rad]

% Resoluci� de les equacions CAS HIPERB�LIC
if isreal(e)==0 || isreal(a)==0 || isreal(theta1)==0 || isempty(e)==1 || isempty(a)==1 || isempty(theta1)==1 || e>1 || e<-1
    [e,a,theta1] = hyperbolic_orbit(r1,r2,deltat,deltatheta);
end

% Correcci� per si surt e negatiu
if e<0
    e = -e;
    theta1 = theta1+pi;
end
theta1 = wrapTo2Pi(theta1);

% C�lcul de la inclinaci� per trigonometria esf�rica
A = asin(cos(beta2)*sin(deltalambda)/sin(deltatheta)); % [rad]
A = checkTangent(A,cos(beta2)*sin(deltalambda),sin(deltatheta));
i = acos(sin(A)*cos(beta1)); % [rad]
l = asin(tan(beta1)/tan(i)); % [rad]

if i>pi/2
    i = i-pi;
end

if beta1==0 && beta2<0 && i>0
    i = -i;
end

if beta1<0 && l>0
    l = -l;
end

sigma = atan(tan(beta1)/cos(A)); % [rad]
sigma = checkTangent(sigma,tan(beta1),cos(A));
Omega = lambda1-l; % [rad]
w = 2*pi-(theta1-sigma); % [rad]

if(i<0)
    i = abs(i);
    Omega = Omega+pi;
    w = w+pi;
end

theta1 = rad2deg(theta1); % [�]
i = rad2deg(i); % [�]
Omega = rad2deg(wrapTo2Pi(Omega)); % [�]
w = rad2deg(wrapTo2Pi(w)); % [�]

end