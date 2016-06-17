close all
clear variables

a0 = 0.4;
da = 0.001;
af = 0.4;
ak = a0:da:af;
T = 20;
% L = pi/sqrt(1-2*a);                                     
% lambda = sqrt(8*a*(1-2*a));                            
% Omega = 2*sqrt(1-2*a);                                  
xj = [0, 0, 0, 0, 0];

for j = 50:50
    for i = 1:length(ak);
        a = ak(i);
        tj = [-1/25*j, 0, 0, 0, 0];
        [PSI, x, t] = calcDarboux(2, a, 256, 1000, T, xj, tj);
        maxima = regions2(PSI, x, t);
    end
end

plot(ak, -time, '-'); grid on;
%figure;
%L = pi./sqrt(1-2*ak);
%plot(L, -time, '.');