function [PSI, xo, to] = calcDarboux(n, a, R, tmax, xs, ts, Nx, Nt, mode, Lk, mult)

if length(a) == 1 && R == 2; % This is the usual default mode
    for k = 1:n
        a(k) = k^2*(a(1)-1/2)+1/2;
    end
elseif length(a) == 1 && R ~= 2;
    if n > 2
        uiwait(warndlg('Ratio not equal to 2 and order > 2. Feature not yet implemented')); 
        return
    end 
    a2 = R^2*(a(1) - 0.5) + 0.5;
    a = [a(1), a2];
end

[~, D] = rat(R);
switch mode
    case 'manual'
        L = Lk;
    case 'periodic'
        L = mult*pi/sqrt(1-2*a(1));                                             % Periodic length
end

l = 1i*sqrt(2*a);
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