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

function [psi] = T1(psi, dt, k2, V, x)
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

pot = V(psi, x);                         % Calculate potential
psi = exp(-1i * pot .* dt).*psi;         % Nonlinear calculation %%% CHANGE 2/2 to dt/2

psi = fftshift(fft(psi));                % FFT
psi = exp(-1i * dt * k2/2).*psi;         % Linear calculation %%%% CHANGE dt*2 to dt
psi = ifft(fftshift(psi));               % Inverse FFT

end