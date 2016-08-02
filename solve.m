function [PSI, x, t] = solve(dt, Nx, Tmax, Lx, mult, V, psi_0, method)
% SOLVE: main function for RogueLab's nonlinear Schrodinger equation solver
% INPUT:
%       dt: temporal step
%       Nx: Number of fourier modes/x-nodes
%       Tmax: maximum simulation time
%       Lx: Box size (usually periodic size for ABs).
%       mult: Box size multiple
%       V: Function handle containing potential
%       psi_0: initial wave function
%       method: Type of algorithm. See below.
% OUTPUT:
%       PSI: spatiotemporal wavefunction matrix
%       x: transverse variable
%       t: temporal variables

hWaitBar = waitbar(0,'Preparing Solver');

Nt = Tmax/dt;                           % Number of temporal nodes
dx = Lx/Nx;                             % Spatial step size
x = (-Nx/2:1:Nx/2-1)'*dx;               % Spatial grid points
t = (0:dt:Tmax).';                      % Temporal grid points
psi = psi_0;                            % Find initial condition
k = 2*(-Nx/2:1:Nx/2-1)'*pi/Lx;          % Wave number
k2 = k.^2;                              % Squares of wavenumbers
PSI = zeros(length(t), length(x));      % Matrix to save whole simulatin
PSI(1, :) = psi;                        % Save first step

%E = energy(psi, k2, Nx, V);            % Initial energy
waitbar(0, hWaitBar, 'Solving: 0%');
for j = 1:Nt                          % Start time evolution
    if strcmp(method, 'T1')             % First order
        psi = T1(psi, dt, k2, V, x); 
    elseif strcmp(method, 'T2')         % Second order
        psi = T2(psi, dt, k2, V, x, mult); 
    elseif strcmp(method, 'T4M')         % Fourth order
        psi = T4M(psi, dt, k2, V, x, mult);
    elseif strcmp(method, 'T4S')      % Fourth order symplectic
        psi = T4S(psi, dt, k2, V, x, mult);
    elseif strcmp(method, 'T6M')         % Sixth order multiproduct
        psi = T6M(psi, dt, k2, V, x, mult);
    elseif strcmp(method, 'T6S')      % Sixth order symplectic
        psi = T6S(psi, dt, k2, V, x, mult);
    elseif strcmp(method, 'T8M')         % Eighth order multiproduct
        psi = T8M(psi, dt, k2, V, x, mult);   
    elseif strcmp(method, 'T8S')      % Eighth order symplectic
        psi = T8S(psi, dt, k2, V, x, mult);
    else                                % Error
        error('Method not recognized. Please consult documentation')
    end
        
    PSI(j+1, :) = psi;                 % Save result of current time step
    
    if ~mod(j/Nt*100, 5)               % Update waitbar every 5% 
        waitbar(j/Nt, hWaitBar, sprintf('Solving: %d%%', j/Nt*100));
    end
end

close(hWaitBar)

end
