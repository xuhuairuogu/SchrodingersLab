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

function [PSI, x, t] = solve(dt, Nx, Tmax, Lx, mult, V, psi_0, method)

hWaitBar = waitbar(0,'Preparing Solver');

Nt = Tmax/dt;                           % Number of temporal nodes
dx = Lx/Nx;                             % Spatial step size
x = (-Nx/2:1:Nx/2-1)'*dx;               % Spatial grid points
t = (0:dt:Tmax).';                          % Temporal grid points
psi = psi_0;                            % Find initial condition
k = 2*(-Nx/2:1:Nx/2-1)'*pi/Lx;          % Wave number
k2 = k.^2;                              % Squares of wavenumbers
PSI = zeros(length(t), length(x));      % Matrix to save whole simulatin
PSI(1, :) = psi;                        % Save first step

%E = energy(psi, k2, Nx, V);            % Initial energy
waitbar(0, hWaitBar, 'Solving: 0%');
for j = 1:1:Nt                          % Start time evolution
    
    if strcmp(method, 'T1')             % First order
        psi = T1(psi, dt, k2, V, x); 
    elseif strcmp(method, 'T2')             % Second order
        psi = T2(psi, dt, k2, V, x, mult); 
        
    elseif strcmp(method, 'T4')         % Fourth order
        psi1 = T2(psi, dt/2, k2, V, x, mult);   
        psi1 = T2(psi1, dt/2, k2, V, x, mult); 
        psi1 = 4/3 * psi1;

        psi2 = T2(psi, dt, k2, V, x, mult); 
        psi2 = -1/3*psi2;

        psi = psi1 + psi2;
        
    elseif strcmp(method, 'T4_NS')      % Fourth order no subtraction
        psi = T4_NS(psi, dt, k2, V, x, mult);
        
    elseif strcmp(method, 'T6')         % Sixth order
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
        
    elseif strcmp(method, 'T6_NS')      % Sixth order no subtraction
        psi = T6_NS(psi, dt, k2, V, x, mult);
        
    elseif strcmp(method, 'T8')         % Eighth order
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
        
    elseif strcmp(method, 'T8_NS')      % Eighth order no subtraction
        psi = T8_NS(psi, dt, k2, V, x, mult);
        
    else                                % Error
        error('Method not recognized. Please consult documentation')
    end
        
    PSI(j+1, :) = psi;
    if ~mod(j/Nt*100, 5)   
        waitbar(j/Nt, hWaitBar, sprintf('Solving: %d%%', j/Nt*100));
    end
end

% Plot results
%waitbar(j/Nt, hWaitBar, 'Preparing Density Plot');
%maxima = regions(PSI, x, t);
% t_c = maxima(2,1);
% disp(t_c);

% maxima = regions(PSI, x, t);
% disp(maxima);
% pk_max = max(max(abs(PSI).^2));
% pk_min = min(min(abs(PSI).^2));
% densityPlot(abs(PSI).^2, x, t, ceil(Nt/1000), ceil(Nx/256), handles.axes2); colormap('jet'); caxis([pk_min pk_max]);
% % ylim([t_c-3, t_c+3]);
% % title('Numerical');
% PSI_k = log(abs(fft(PSI'))/Nx);
% fourierPlot(PSI_k', t, 7, mult, handles.axes3);
% % xlim([t_c-3, t_c+3]);
% % ylim([-10, 2]);
% % title('Numerical');
% %PSI_k = fftshift(log(abs(fft(PSI'))/Nx), 1);
% %densityPlot(PSI_k', k, t, Nt/1000, Nx/256, handles.axes3); colormap('jet'); 
% initialPlot(PSI(1, :), x, handles.axes1);
close(hWaitBar)
% 
% % Some special purpose functions
% % b_plot(PSI, Nx, t, 6, maxima(1,1));
% % recon(Nx, 5000, max(t), maxima(1,1));
% % [~,~,shift] = ab(PSI, x, t, max(max(abs(psi).^2)), 3/8);      
% % energy(PSI, t, k2, Nx, V, dt);
end
