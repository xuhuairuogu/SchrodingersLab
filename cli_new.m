% MAIN CLI SOLVER FOR SchrodingerMAT
clear variables
close all

Nx = 2^7;                   % Number of modes
dt = 0.001;   % Temporal spacing delta t
Tmax = 80;     % Max time for simulation
Nt = Tmax/dt;

V = @(psi, x) (-1*abs(psi).^2);          % Potential
method = 'T2';              % Algorithm to use, this means 4th order algorithm with multi-product integrators

gamma0 = 20; % Other parameters for Absorbing BC
alpha = 1;   % Other parameter for absorbing BC
gamma = @(x, gamma_0, alpha, L) gamma_0*(sech(alpha*(x-L/2)).^2 + sech(alpha*(x+L/2)).^2);        % Absorbing boundary conditions function
absorption = struct('useAbsorbingBC', 0, 'gamma', gamma, 'gamma0', gamma0, 'alpha', alpha);       %absorbing BC structure

f1 = figure;
h1 = axes();
f2 = figure;
h2 = axes();
f3 = figure;
h3 = axes();
handles = struct('axes1', h1, 'axes2', h2, 'axes3', h3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USUAL INITIAL COND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A1 = exp(-6);
% A0 = sqrt(1-2*A1^2);
% a = 0.386851;
% Omega = 2*sqrt(1-2*a);
% L = mult*pi/sqrt(1-2*a);
% dx = L/Nx;                             % Spatial step size
% x = (-Nx/2:1:Nx/2-1)'*dx;   
% psi_0 = A0 + 2*A1*cos(Omega*x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Darboux Transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%psi_0 = calcDarboux(nu, Nx, Nt, 6.3*2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Triplet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L = 20;
% beta = 0;
% gamma = 0;
% t0 = 0;
% psi_0 = triplet(beta, gamma, Nx, 200, L, t0, t0+2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARTIFICIALLY SMOOTH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L = 100;
% M = 1;
% dx = L/Nx;                             % Spatial step size
% xo = (-Nx/2:1:Nx/2-1)'*dx;               % Spatial grid points
% psi_i = -0.43 + 0.44*sech(xo/M);
% psi_r = -0.84-0.19*sech(xo/M).^2;
% psi_0 = psi_r' + 1i*psi_i';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARTIFICAL 2 AB's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%a1 = 0.375;
% nu = 2;
% a1= 0.5*(1-1/nu^2);
% L = 5*pi/sqrt(1-2*a1);
% v1 = 0;
% psi_0 = analytical_ab(v1, a1, Nx, 100, L, 0, 5);
% psi_02 = analytical_ab(-v1, a1, Nx/2, 100, L/2, 0, 5);
% psi_0 = [psi_01, psi_02];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_0 = 0.380;
da = 0.0001;
a_f = 0.385;

ticM
for a =a_0:da:a_f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYTICAL AB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nu = 2;
% a = 0.5*(1-1/nu^2);
% a = 0.383;
Omega = 2*sqrt(1-2*a);
mult = 1;     % KEEP AT ONE EXCEPT WITH NL TALBOT CARPETS IN WHICH CASE UNCOMMENT FROM T2
L = mult*pi/sqrt(1-2*a);
dx = L/Nx;     
x = (-Nx/2:1:Nx/2-1)'*dx;  
t_s = 0;
x_s = x;
lambda = Omega*sqrt(1-Omega^2/4);
psi_s = ((1-4*a)*cosh(lambda*t_s)+sqrt(2*a)*cos(Omega*x_s)+1i*lambda*sinh(lambda*t_s))./(sqrt(2*a)*cos(Omega*x_s)-cosh(lambda*t_s));
psi_0 = psi_s;

solve(dt, Nx, Tmax, L, mult, V, psi_0, absorption, method, handles);
title(h2, sprintf('a = %f', a));
title(h3, sprintf('a = %f', a));
print(f2, '-dpng2', sprintf('test/intensity_a_%f.png', a))
print(f3, '-dpng2', sprintf('test/spectrum_a_%f.png', a))
disp(['Current a = ', num2str(a)]);
hold off;
end
toc;
load handel;
p = audioplayer(y, Fs);
play(p, [1 (get(p, 'SampleRate') * 3)]);