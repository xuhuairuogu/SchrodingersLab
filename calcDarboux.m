function psi_0 = calcDarboux(nu, Nx, Nt, tmax)

n = 1;                                                          % Order of rogue wave
p = 1;                                                          % For enumeration. Leave as 1

%nu = 3;                                                         % For first order waves, number of unstable modes + 1

%a = 0.5*(1-1/(nu*nu));                                          % Parameter a

%l = 1i*sqrt(2*a);                                               % Relating eigenvalue to parameter a
a = 0.375;
l = [-0.9 + 1i*sqrt(2*a), 0.9 + 1i*sqrt(2*a)];
L = nu*pi/sqrt(1-2*a);                                             % Periodic length
kappa = 2*sqrt(1+l.^2);                                         % Principal wave number
ts = [0, 0];                                                    % t-shifts
xs = [0, 0];                                                    % x-shifts

psi = cell(n, 1);                                               % container for all order waves
dt = tmax/Nt;
to = (-Nt/2:1:Nt/2-1)'*dt;                                     % propagation variable
dx = L/Nx;                             % Spatial step size
xo = (-Nx/2:1:Nx/2-1)'*dx;   
[x, t] = meshgrid(xo, to);                                      % Mesh for calculation

[~, ~, psi] = calc_rs(n, p, x, t, xs, ts, l, kappa, psi);       % Start recursion
figure;
h1 = axes();
PSI = psi{n};
densityPlot(abs(PSI).^2, xo, to, Nt/2000, Nx/512, h1); colormap('jet');        % Plot result

psi_0 = PSI(1, :);
maxima = regions(PSI, xo, to);
disp(maxima);

%figure;
%plot(x, abs(psi_0).^2);

%max(abs(psi_0).^2)

end