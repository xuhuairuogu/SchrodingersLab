% MAIN CLI SOLVER FOR SchrodingerMAT
clear variables
close all

nu = 4;
q = 0.5*(1-1/(nu*nu));
Omega = 2*sqrt(1-2*q);

A1 = 1e-2;
A0 = sqrt(1-2*A1^2);
Nx = 256;                   % Number of modes
%k = 1.2;
%L = 4*4/sqrt(2*k)*ellipk(sqrt(k-1)/sqrt(2*k));                      % Box size
%L = 160;
L = pi/sqrt(1-2*q);

dt = 0.001;   % Temporal spacing delta t
Tmax = 40;     % Max time for simulation
V = @(psi, x) (-1*abs(psi).^2);          % Potential
gamma = @(x, gamma_0, alpha, L) gamma_0*(sech(alpha*(x-L/2)).^2 + sech(alpha*(x+L/2)).^2);        % Absorbing boundary conditions function
method = 'T2';              % Algorithm to use, this means 4th order algorithm with multi-product integrators

dx = L/Nx;                             % Spatial step size
x = (-Nx/2:1:Nx/2-1)'*dx;           
t = 0;

psi_0 = A0 + 2*A1*cos(Omega*x)';
%psi_0 = 1- 4./(1+4*x.^2)';
% psi_0 = AB_je(x, 0, k)';
% a1 = 4/9;
% Omega1 = 2*sqrt(1-2*a1);
% chi1 = 1/2*acos(Omega1/2);
% lambda1 = sqrt(8*a1*(1-2*a1));
% 
% a2 = 1/2-(2*Omega1)^2/8;
% Omega2 = 2*sqrt(1-2*a2);
% chi2 = 1/2*acos(Omega2/2);
% lambda2 = sqrt(8*a2*(1-2*a2));
% 
% f2 = Omega2*x;
% gamma_r2 = 1/2*f2 + chi2;
% gamma_s2 = 1/2*f2 - chi2;
% 
% r2 = -1i*sqrt(2)*(cos(gamma_r2) - sin(gamma_r2));
% s2 =     sqrt(2)*(cos(gamma_s2) + sin(gamma_s2));
% 
% B = sin(Omega1*x)*lambda1/2;
% A = sqrt(2*a1)*cos(Omega1*x) - 2*a1;
% C = cos(Omega1*x)*sqrt(2*a1) - 1;
% A = A./C;
% B = B./C;
% 
% s12 =  1i*sqrt(2*a2)*s2 -  (A.*r2 - 1i*B.*s2);
% r12 =  1i*sqrt(2*a2)*r2  + (A.*s2 - 1i*B.*r2);
% 
% psi_1 = ((1-4*a1)+sqrt(2*a1)*cos(Omega1*x))./(sqrt(2*a1)*cos(Omega1*x)-1);
% psi_12 = -psi_1 + 1*(-4*1i*sqrt(2*a2)*s12.*conj(r12))./(abs(r12).^2 + abs(s12).^2);

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
solve(dt, Nx, Tmax, L, V, psi_0', absorption, method, handles);