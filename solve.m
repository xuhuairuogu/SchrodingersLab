% Copyright 2015, Omar Ashour.
% This sourcecode is available from <https://github.com/oashour/HighNLSE/>
%
% This file is part of HighNLSE.
% 
% HighNLSE is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% HighNLSE is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with HighNLSE.  If not, see <http://www.gnu.org/licenses/>.

function solve(dt, Nx, Tmax, Lx, V, psi_0, method, handles)

hWaitBar = waitbar(0,'Preparing Solver');

Nt = Tmax/dt;                           % Number of temporal nodes
dx = Lx/Nx;                             % Spatial step size
x = (-Nx/2:1:Nx/2-1)'*dx;               % Spatial grid points
t = 0:dt:Tmax;                          % Temporal grid points
psi = eval(psi_0);                      % Find initial condition
k = 2*(-Nx/2:1:Nx/2-1)'*pi/Lx;          % Wave number
k2 = k.^2;                              % Squares of wavenumbers
PSI = zeros(length(t), length(x));      % Matrix to save whole simulatin
PSI(1, :) = psi;                        % Save first step

%E = energy(psi, k2, Nx, V);            % Initial energy
waitbar(0, hWaitBar, 'Solving: 0%%');
for j = 1:1:Nt                          % Start time evolution
    
    if strcmp(method, 'T2')             % Second order
        psi = T2(psi, dt, k2, V, x); 
        
    elseif strcmp(method, 'T4')         % Fourth order
        psi1 = T2(psi, dt/2, k2, V, x);   
        psi1 = T2(psi1, dt/2, k2, V); 
        psi1 = 4/3 * psi1;

        psi2 = T2(psi, dt, k2, V); 
        psi2 = -1/3*psi2;

        psi = psi1 + psi2;
        
    elseif strcmp(method, 'T4_NS')      % Fourth order no subtraction
        psi = T4_NS(psi, dt, k2, V);
        
    elseif strcmp(method, 'T6')         % Sixth order
        psi1 = T2(psi, dt/3, k2, V);   
        psi1 = T2(psi1, dt/3, k2, V); 
        psi1 = T2(psi1, dt/3, k2, V);
        psi1 = 81/40 * psi1;

        psi2 = T2(psi, dt/2, k2, V); 
        psi2 = T2(psi2, dt/2, k2, V);
        psi2 = -16/15*psi2;

        psi3 = T2(psi, dt, k2, V); 
        psi3 = 1/24*psi3;

        psi = psi1 + psi2 + psi3;
        
    elseif strcmp(method, 'T6_NS')      % Sixth order no subtraction
        psi = T6_NS(psi, dt, k2, V);
        
    elseif strcmp(method, 'T8')         % Eighth order
        psi1 = T2(psi, dt/4, k2, V);   
        psi1 = T2(psi1, dt/4, k2, V); 
        psi1 = T2(psi1, dt/4, k2, V);
        psi1 = T2(psi1, dt/4, k2, V);
        psi1 = 1024/315 * psi1;

        psi2 = T2(psi, dt/3, k2, V); 
        psi2 = T2(psi2, dt/3, k2, V);
        psi2 = T2(psi2, dt/3, k2, V);
        psi2 = -729/280*psi2;

        psi3 = T2(psi, dt/2, k2, V); 
        psi3 = T2(psi3, dt/2, k2, V);
        psi3 = 16/45*psi3;

        psi4 = T2(psi, dt, k2, V); 
        psi4 = -1/360*psi4;

        psi = psi1 + psi2 + psi3 + psi4;
        
    elseif strcmp(method, 'T8_NS')      % Eighth order no subtraction
        psi = T8_NS(psi, dt, k2, V);
        
    else                                % Error
        error('Method not recognized. Please consult documentation')
    end
        
    PSI(j+1, :) = psi;
    if ~mod(j/Nt*100, 5)   
        waitbar(j/Nt, hWaitBar, sprintf('Solving: %d%%', j/Nt*100));
    end
end

% Plot results
waitbar(j/Nt, hWaitBar, 'Preparing Density Plot');
densityPlot(PSI, x, t, dt, dx, 200, handles.axes2);    
colormap('jet');
PSI_k = 20*log10(abs(fft(PSI'))/Nx/Nx);
max(max(PSI_k))
min(min(PSI_k))
densityPlot(fftshift(PSI_k,1)', k, t, dt, dx, 200, handles.axes3);  
colormap('jet');
initialPlot(PSI(1, :), x, handles.axes1);
close(hWaitBar)

% Some special purpose functions
maxima = regions(PSI, x, t);
% b_plot(PSI, Nx, t, 6, maxima(1,1));
% recon(Nx, 5000, max(t), maxima(1,1));
% [~,~,shift] = ab(PSI, x, t, max(max(abs(psi).^2)), 3/8);      
% energy(PSI, t, k2, Nx, V, dt);
end
