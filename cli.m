% MAIN CLI SOLVER FOR SchrodingerMAT
clear variables
close all

L = 80;                      % Box size
Nx = 1024;                   % Number of modes

psi = @(x) (1 - 4./(1+4*x.^2));         % Initial conditions
dt = 0.001;   % Temporal spacing delta t
Tmax = 20;     % Max time for simulation
V = @(psi, x) (-1*abs(psi).^2);          % Potential
gamma = @(x, gamma_0, alpha, L) gamma_0*(sech(alpha*(x-L/2)).^2 + sech(alpha*(x+L/2)).^2);        % Absorbing boundary conditions function
method = 'T4';              % Algorithm to use, this means 4th order algorithm with multi-product integrators

% Prepare figures to display results
figure;
h1 = axes();
figure;
h2 = axes();
figure;
h3 = axes();
handles = struct('axes1', h1, 'axes2', h2, 'axes3', h3);

gamma0 = 20; % Other parameters for Absorbing BC
alpha = 1;   % Other parameter for absorbing BC

absorption = struct('useAbsorbingBC', 0, 'gamma', gamma, 'gamma0', gamma0, 'alpha', alpha);       %absorbing BC structure
solve(dt, Nx, Tmax, L, V, psi, absorption, method, handles); 
