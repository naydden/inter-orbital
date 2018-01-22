function [e,a,theta] = hyperbolic_orbit(r1,r2,deltat,deltatheta)

resta = 1000;
d = 1e-6;
theta = deg2rad(-60); % valor inicial que es dóna a theta
eant = 1000;

while(resta>d)
    
    e = (norm(r2)-norm(r1))/(norm(r1)*cos(theta)-norm(r2)*cos(theta+deltatheta));
    a = norm(r1)*(1+e*cos(theta))/(e^2-1);
    delta = 365.25*a^(3/2)*(e*sqrt(e^2-1)*sin(theta+deltatheta)/(1+e*cos(theta+deltatheta))-log(abs((tan((theta+deltatheta)/2)+sqrt((e+1)/(e-1)))/(tan((theta+deltatheta)/2)-sqrt((e+1)/(e-1)))))-e*sqrt(e^2-1)*sin(theta)/(1+e*cos(theta))+log(abs((tan(theta/2)+sqrt((e+1)/(e-1)))/(tan(theta/2)-sqrt((e+1)/(e-1))))))/(2*pi);
    
    resta = abs(e-eant);
    eant = e;
    
    if resta>d
        theta = theta+(deltat-delta)/1000;
    end
    
end
end