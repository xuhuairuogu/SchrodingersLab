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

function [error, E, ke, pe] = energy(psi, t, k2, Nx, gamma, dt)
% ENERGY: finds the change in energy dE throughout the whole simulation
% INPUT:
%       psi: a matrix with the wave function evolution in space and time.
%            rows are time, columns are space.
%       t: time array
%       k2: square of wave number
%       Nx: number of Fourier modes/spatial nodes
%       gamma: strength of the nonlinearity
%       dt: temporal separation
% OUTPUT:
%       error: integrated energy error
%       E: array with energy at every time step
%       ke: array with kinetic energy at every time step
%       pe: array with potential energy at every time step


pe = zeros(1, length(t));                           % Preallocate pe
ke = zeros(1, length(t));                           % Preallocate ke

for i = 1:length(t)                                 % loop over time
    psi2 = abs(psi(i, :)').*abs(psi(i, :)');        % square psi
    pe(i) = sum(0.5*gamma*psi2.*psi2)/sum(psi2);    % calculate PE

    psi_k = fftshift(fft(psi(i, :)'))/Nx;           % Calculate normalized FFT
    psi_k2 = abs(psi_k).*abs(psi_k);                % Square normalized FFT
    ke(i) =  sum(0.5*k2.*psi_k2)/sum(psi_k2);       % Calculate ke
end

E = ke + pe;                                        % Sum PE and KE to get E

E0 = E(1);                                          % Initial energy
dE = E-E0;                                          % Energy error

figure
plot(t, dE, '-b', 'LineWidth', 2)                   % Plot
xlabel('t'); ylabel('E-E0');                        % Legend
title(sprintf('dt = %0.4f', dt))                    % Title

error = sum(dE)*dt;                                 % integrated energy error
end

