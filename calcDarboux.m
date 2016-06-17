function [PSI, xo, to] = calcDarboux(n, a, Nx, Nt, tmax, xs, ts)

%a = 0.5*(1-1/(nu*nu));                                          % Parameter a

a1 = a;
a2 = 4*a-3/2;
a3 = 9*a - 4;
a4 = 16*a - 15/2;
a5 = 25*a - 12;
l = [1i*sqrt(2*a1), 1i*sqrt(2*a2), 1i*sqrt(2*a3), 1i*sqrt(2*a4), 1i*sqrt(2*a5)];
L = pi/sqrt(1-2*a);                                             % Periodic length
kappa = 2*sqrt(1+l.^2);                                         % Principal wave number

psi = cell(n, 1);                                               % container for all order waves
if Nt > 0
    dt = tmax/Nt;
    to = (-Nt/2:1:Nt/2-1)'*2*dt;                                     % propagation variable
else 
    to = tmax;
end

if Nx > 0
    dx = L/Nx;                             % Spatial step size
    xo = (-Nx/2:1:Nx/2-1)'*dx; 
else
    xo = 0;
end
[x, t] = meshgrid(xo, to);                                      % Mesh for calculation

[~, ~, psi] = calc_rs(n, 1, x, t, xs, ts, l, kappa, psi);       % Start recursion
PSI = psi{n};
% maxima = regions(PSI, xo, to);
% disp(maxima);
% %max(abs(psi_0).^2)
end