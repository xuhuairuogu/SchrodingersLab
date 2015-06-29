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

function [psi] = T2(psi, dt, k2, gamma)
% T2:
%   This function calculates one time step using a second order split step
%   algorithm. 
%   OUTPUT:
%       psi: evolution after one time step.
%   INPUT:
%       psi: function in spatial domain (not fft'd).
%       dt:  temporal separation between nodes.
%       k2:  square of wave numbers. Do not fft shift this.
%       gamma: nonlinearity.
%
%   FORTRAN CODE:
%       epot=agp*cdabs(psi)**2
%       phi = cdexp( -im * 0.5d0*dt * epot )*psi
%       call fft (phi, m, 0 )
%       phi =cdexp( -im * dt * k2/2 )*phi
%       call fft (phi, m, 1 )
%       epot=agp*cdabs(phi)**2
%       psi = cdexp( -im * 0.5d0*dt * epot )*phi

psi2 = gamma*abs(psi).^2;                % Calculate potential
psi = exp(-1i * dt/2 * psi2).*psi;       % Nonlinear calculation
psi = fftshift(fft(psi));                % FFT
psi = exp(-1i * dt * k2/2).*psi;         % Linear calculation
psi = ifft(fftshift(psi));               % Inverse FFT
psi2 = gamma*abs(psi).^2;                % Calculate potential
psi = exp(-1i * dt/2 * psi2).*psi;       % nonlinear calculation

end

