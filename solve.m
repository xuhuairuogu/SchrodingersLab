function [PSI, x, t, k] = solve(dt, Nx, Tmax, Lx, mult, V, psi_0, algorithm)
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

algList = {'T1', 'T2', 'T4M', 'T4S', 'T6M', 'T6S', 'T8M', 'T8S'};
% Select Algorithm
if ~any(strcmp(algorithm,algList))
    error('Algorithm unrecognized. Please check documentation')
else
    Tn = str2func(algorithm);
end

waitbar(0, hWaitBar, 'Solving: 0%');
for j = 1:Nt                          % Start time evolution
    PSI(j+1, :) = Tn(PSI(j,:).', dt, k2, V, x, mult);  % Evolve 1 step
    
    if ~mod(j/Nt*100, 5)               % Update waitbar every 5% 
        waitbar(j/Nt, hWaitBar, sprintf('Solving: %d%%', j/Nt*100));
    end
end

close(hWaitBar)

function [psi] = T1(psi, dt, k2, V, x, ~)
% T1:
%   This function calculates one time step using a first order split step
%   algorithm. 
%   OUTPUT:
%       psi: evolution after one time step.
%   INPUT:
%       psi: function in spatial domain (not fft'd).
%       dt:  temporal separation between nodes.
%       k2:  square of wave numbers. Do not fft shift this.
%       V:   potential

pot = V(psi, x);                         % Calculate potential
psi = exp(-1i * pot .* dt).*psi;         % Nonlinear calculation

psi = fftshift(fft(psi));                % FFT
psi = exp(-1i * dt * k2/2).*psi;         % Linear calculation 
psi = ifft(fftshift(psi));               % Inverse FFT

function [psi] = T2(psi, dt, k2, V, x, mult)
% T2:
%   This function calculates one time step using a second order split step
%   algorithm. 
%   OUTPUT:
%       psi: evolution after one time step.
%   INPUT:
%       psi: function in spatial domain (not fft'd).
%       dt:  temporal separation between nodes.
%       k2:  square of wave numbers. Do not fft shift this.
%       V:   potential
Nx = length(x);

pot = V(psi, x);                        % Calculate potential
psi = exp(-1i * dt/2 * pot).*psi;        % Nonlinear calculation

psi = fftshift(fft(psi));                % FFT
psi = exp(-1i * dt * k2/2).*psi;         % Linear calculation
psi = ifft(fftshift(psi));               % Inverse FFT

pot = V(psi, x);                         % Calculate potential
psi = exp(-1i * dt/2 * pot).*psi;       % Nonlinear calculation

% See documentation for explanation of what mult does.
% You will almost always want to set mult to 1 unless you are investigating
% Nonlinear talbot carpets formed by AB's. For example, let's say your mult
% is set to 3, i.e. box size is 3*periodic length. This will result in the
% non-triplet modes growing and ruining your carpet, so this sets them to 0
% and kills their growth artificially. Pretty much useful for nothing but
% generating ugly dress patterns for your grandmother.
if mult > 1                             
    psi = fft(psi);              
    for i = 2:Nx/2+1
        if(mod(i-1, mult) ~= 0)
            psi(i) = 0;
            psi(Nx - i  + 2) = 0;
        end
    end
    
    psi = ifft(psi);
end

function psi = T4M(psi, dt, k2, V, x, mult) %#ok<*DEFNU>
% T4M:
%   This function calculates one time step using an fourth order split step
%   multi-product algorithm. 
%   OUTPUT:
%       psi: evolution after one time step.
%   INPUT:
%       psi: function in spatial domain (not fft'd).
%       dt:  temporal separation between nodes.
%       k2:  square of wave numbers. Do not fft shift this.
%       V:   potential (function handle)

psi1 = T2(psi, dt/3, k2, V, x, mult);    
psi1 = T2(psi1, dt/3, k2, V, x, mult);  
psi1 = T2(psi1, dt/3, k2, V, x, mult); 
psi1 = 81/40 * psi1;

psi2 = T2(psi, dt/2, k2, V, x, mult); 
psi2 = T2(psi2, dt/2, k2, V, x, mult); 
psi2 = -16/15*psi2;

psi3 = T2(psi, dt, k2, V, x, mult); 
psi3 = 1/24*psi3;

psi = psi1 + psi2 + psi3;

function psi = T4S(psi, dt, k2, V, x, mult)
% T4S:
%   This function calculates one time step using a fourth order split step
%   symplectic algorithm. 
%   OUTPUT:
%       psi: evolution after one time step.
%   INPUT:
%       psi: function in spatial domain (not fft'd).
%       dt:  temporal separation between nodes.
%       k2:  square of wave numbers. Do not fft shift this.
%       V:   potential
s = 2^(1/3);
os = 1/(2-s);

ft = os;                % Forward factor
bt = -s*os;             % Backward factor
 
psi = T2(psi, ft*dt, k2, V, x, mult);
psi = T2(psi, bt*dt, k2, V, x, mult);
psi = T2(psi, ft*dt, k2, V, x, mult);

function psi = T6M(psi, dt, k2, V, x, mult)
% T6M:
%   This function calculates one time step using an sixth order split step
%   multi-product algorithm. 
%   OUTPUT:
%       psi: evolution after one time step.
%   INPUT:
%       psi: function in spatial domain (not fft'd).
%       dt:  temporal separation between nodes.
%       k2:  square of wave numbers. Do not fft shift this.
%       V:   potential (function handle)

psi1 = T2(psi, dt/3, k2, V, x, mult);    
psi1 = T2(psi1, dt/3, k2, V, x, mult);  
psi1 = T2(psi1, dt/3, k2, V, x, mult); 
psi1 = 81/40*psi1;

psi2 = T2(psi, dt/2, k2, V, x, mult); 
psi2 = T2(psi2, dt/2, k2, V, x, mult); 
psi2 = -16/15*psi2;

psi3 = T2(psi, dt, k2, V, x, mult); 
psi3 = 1/24*psi3;

psi = psi1 + psi2 + psi3;

function psi = T6S(psi, dt, k2, V, x, mult)
% T6S:
%   This function calculates one time step using a sixth order split step
%   symplectic algorithm. 
%   OUTPUT:
%       psi: evolution after one time step.
%   INPUT:
%       psi: function in spatial domain (not fft'd).
%       dt:  temporal separation between nodes.
%       k2:  square of wave numbers. Do not fft shift this.
%       V:   potential
s = 2^(1/5);
os = 1/(2-s);

ft = os;                % Forward Step
bt = -s*os;             % Backward Step

psi = T4S(psi, ft*dt, k2, V, x, mult);
psi = T4S(psi, bt*dt, k2, V, x, mult);
psi = T4S(psi, ft*dt, k2, V, x, mult);

function psi = T8M(psi, dt, k2, V, x, mult)
% T8M:
%   This function calculates one time step using an eighth order split step
%   multi-product algorithm. 
%   OUTPUT:
%       psi: evolution after one time step.
%   INPUT:
%       psi: function in spatial domain (not fft'd).
%       dt:  temporal separation between nodes.
%       k2:  square of wave numbers. Do not fft shift this.
%       V:   potential

psi1 = T2(psi, dt/4, k2, V, x, mult);   
psi1 = T2(psi1, dt/4, k2, V, x, mult);
psi1 = T2(psi1, dt/4, k2, V, x, mult);
psi1 = T2(psi1, dt/4, k2, V, x, mult);
psi1 = 1024/315 * psi1;

psi2 = T2(psi, dt/3, k2, V, x, mult); 
psi2 = T2(psi2, dt/3, k2, V, x, mult);
psi2 = T2(psi2, dt/3, k2, V, x, mult);
psi2 = -729/280*psi2;

psi3 = T2(psi, dt/2, k2, V, x, mult); 
psi3 = T2(psi3, dt/2, k2, V, x, mult);
psi3 = 16/45*psi3;

psi4 = T2(psi, dt, k2, V, x, mult);
psi4 = -1/360*psi4;

psi = psi1 + psi2 + psi3 + psi4;

function [psi] = T8S(psi, dt, k2, V, x, mult)
% T8S:
%   This function calculates one time step using an eighth order split step
%   symplectic algorithm. 
%   OUTPUT:
%       psi: evolution after one time step.
%   INPUT:
%       psi: function in spatial domain (not fft'd).
%       dt:  temporal separation between nodes.
%       k2:  square of wave numbers. Do not fft shift this.
%       V:   potential

s = 2^(1/7);
os = 1/(2-s);

ft = os;            % Forward Step
bt = -s*os;         % Backward Step

psi = T6S(psi, ft*dt, k2, V, x, mult);
psi = T6S(psi, bt*dt, k2, V, x, mult);
psi = T6S(psi, ft*dt, k2, V, x, mult);