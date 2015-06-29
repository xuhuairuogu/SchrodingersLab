% High NLSE, a MATLAB solver for the nonlinear Schrodinger Equation.
% Coded by Omar Ashour, Texas A&M University at Qatar.
% Last revision: 29/6/2015.

close all
clear all

dt = 0.01;                              % Temporal step size
Nx = 2^7;                               % Number of fourier modes/spatial nodes
Tmax = 50;                              % Maximum time to run simulation
Nt = Tmax/dt;                           % Number of temporal nodes
intr = 10;                              % Interval for caputring image for graph
Lx = 2*pi;                              % box size = [-Lx/2, Lx/2)
psi_0 = '1 + 1e-16*cos(3*x)';           % Initial wave function
gamma = -1;                             % Strength of nonlinearity
method = 'T2';

dx = Lx/Nx;                             % Spatial step size
x = (-Nx/2:1:Nx/2-1)'*dx;               % Spatial grid points
t = 0:dt:Tmax;                          % Temporal grid points
psi = eval(psi_0);                      % Find initial condition
k = 2*(-Nx/2:1:Nx/2-1)'*pi/Lx;          % Wave number
k2 = k.^2;                              % Squares of wavenumbers
PSI = zeros(length(t), length(x));      % Matrix to save whole simulatin
PSI(1, :) = psi;                        % Save first step

%E = energy(psi, k2, Nx, gamma);        % Initial energy energy
for j = 1:1:Nt                          % Start time evolution
    
    if strcmp(method, 'T2')             % Second order
        psi = T2(psi, dt, k2, gamma); 
        
    elseif strcmp(method, 'T4')         % Fourth order
        psi1 = T2(psi, dt/2, k2, gamma);   
        psi1 = T2(psi1, dt/2, k2, gamma); 
        psi1 = 4/3 * psi1;

        psi2 = T2(psi, dt, k2, gamma); 
        psi2 = -1/3*psi2;

        psi = psi1 + psi2;
        
    elseif strcmp(method, 'T4_NS')      % Fourth order no subtraction
        psi = T4_NS(psi, dt, k2, gamma);
        
    elseif strcmp(method, 'T6')         % Sixth order
        psi1 = T2(psi, dt/3, k2, gamma);   
        psi1 = T2(psi1, dt/3, k2, gamma); 
        psi1 = T2(psi1, dt/3, k2, gamma);
        psi1 = 81/40 * psi1;

        psi2 = T2(psi, dt/2, k2, gamma); 
        psi2 = T2(psi2, dt/2, k2, gamma);
        psi2 = -16/15*psi2;

        psi3 = T2(psi, dt, k2, gamma); 
        psi3 = 1/24*psi3;

        psi = psi1 + psi2 + psi3;
        
    elseif strcmp(method, 'T6_NS')      % Sixth order no subtraction
        psi = T6_NS(psi, dt, k2, gamma);
        
    elseif strcmp(method, 'T8')         % Eighth order
        psi1 = T2(psi, dt/4, k2, gamma);   
        psi1 = T2(psi1, dt/4, k2, gamma); 
        psi1 = T2(psi1, dt/4, k2, gamma);
        psi1 = T2(psi1, dt/4, k2, gamma);
        psi1 = 1024/315 * psi1;

        psi2 = T2(psi, dt/3, k2, gamma); 
        psi2 = T2(psi2, dt/3, k2, gamma);
        psi2 = T2(psi2, dt/3, k2, gamma);
        psi2 = -729/280*psi2;

        psi3 = T2(psi, dt/2, k2, gamma); 
        psi3 = T2(psi3, dt/2, k2, gamma);
        psi3 = 16/45*psi3;

        psi4 = T2(psi, dt, k2, gamma); 
        psi4 = -1/360*psi4;

        psi = psi1 + psi2 + psi3 + psi4;
        
    elseif strcmp(method, 'T8_NS')      % Eighth order no subtraction
        psi = T8_NS(psi, dt, k2, gamma);
        
    else                                % Error
        error('Method not recognized. Please consult documentation')
    end
    
    disp([j Nt])
        
    PSI(j+1, :) = psi;
end

% Plot results
surf(x,t,abs(PSI).^2, 'EdgeColor', 'none');
ylim([0, max(t)])
xmin = x(1); xmax = x(Nx);
xlim([xmin, xmax])
colorbar('eastoutside')
ylabel('t'); xlabel('x'); zlabel('|\psi|^2');

% Some special purpose functions
% a_plot(PSI, Nx, t, 5);
% recon(37.883, 0.5333, 1.86, Nx, 5000, max(t));
% ab(PSI, x, t, max(max(abs(psi).^2)), 3/8);      
% energy(PSI, t, k2, Nx, gamma, dt);