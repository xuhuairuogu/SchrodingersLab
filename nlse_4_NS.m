function [error, minimum, maximum] = nlse_4_NS(dt)
% nlse_2(dt)
%   Solves the nonlinear schrodinger equation using a second order split.
%   step algorithm.
%   OUTPUT:
%       error: the integrated error in the energy.
%   INPUT:
%       dt: temporal separation between nodes.
%       Nx: Number of fourier modes/spatial nodes
%       Tmax: Maximum time to run simulation
%       Nt: Number of temporal nodes
%       interval: Interval for caputring image for graph
%       Lx: box size = [-Lx/2, Lx/2)
%       psi_0: Initial wave function (string)
%       gamma:  Strength of nonlinearity

Nx = 2^7;                               % Number of fourier modes/spatial nodes
Tmax = 50;                            % Maximum time to run simulation
Nt = Tmax/dt;                           % Number of temporal nodes
interval = 10;                          % Interval for caputring image for graph
Lx = 2*pi;                              % box size = [-Lx/2, Lx/2)
psi_0 = '1 + 1e-6*cos(3*x)';           % Initial wave function
gamma = -1;                             % Strength of nonlinearity

%dx = Lx/Nx;                             % Spatial step size
x = (-Nx/2:1:Nx/2-1)'*Lx/Nx;               % Grid points
t = 0:dt:Tmax;
psi = eval(psi_0); PSI = psi;           % Find initial condition
FULL = psi;
k = 2*(-Nx/2:1:Nx/2-1)'*pi/Lx;          % Wavenumbers
k2 = k.^2;                              % Squares of wave numbers

E = energy(psi, k2, Nx, gamma);         % Initial energy energy


for m = 1:1:Nt                          % Start time evolution

   % T4(dt) = (4/3) T2(dt/2)T(dt/2) ? (1/3) T2(dt)

    psi = T4_NS(psi, dt, k2, gamma);

    E = [E energy(psi, k2, Nx, gamma)];% Energy at current time step
   
    FULL = [FULL psi];
    if rem(m,interval) == 0                      
        PSI = [PSI psi];                % Save results
        disp([m Nt])                    % Display progress
    end                                  
end

% Plot results
figure
d = (0:interval:Nt)';
t2 = d*dt;
surf(x,t2,abs(PSI).^2', 'EdgeColor', 'none')
ylim([0, max(t2)])
xmin = x(1); xmax = x(Nx);
xlim([xmin, xmax])
colorbar('eastoutside')
ylabel('t')
xlabel('x')
zlabel('|\psi|^2')

% Plot energy error
E0 = E(1);
dE = E-E0;
figure
plot(linspace(0, Tmax, Nt+1), dE)
xlabel('Time'); ylabel('E-E0'); title(sprintf('dt = %0.4f', dt))

error = sum(dE)*dt;                     % Find integrated energy error
minimum = min(dE);
maximum = max(dE);

[m1 i1] = max(abs(FULL).^2);
[m11 j1] = max(m1);
k1 = i1(j1);
ab(FULL, x, t, k1, j1, m11);
