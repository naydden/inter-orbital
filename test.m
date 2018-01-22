ro = [8000000: 10000 : 11000000];
load('DATA.mat');
mu = earth.mu;
for i=1:301
    Vo = sqrt(mu/ro(i));
    Vinf(i) = sqrt(3.6054e+03^2+2*Vo^2)-Vo;
end
plot(ro,Vinf)