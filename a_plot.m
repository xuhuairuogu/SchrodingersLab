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

function [] = a_plot(PSI, Nx, t, num)
% FUNCTION: Plots analytical function of r, c, d, e vs actual data
% INPUT:
%       PSI: full spatiotemporal wave function
%       Nx: Number of Fourier modes/spatial nodes
%       t: time array
%       num: number of elements to plot. 5 is recommended.
% See the file recon for more info on the analytical fit.

% Prepare the actual data
J = 10;                                     % Scaling down of number of points
PSI_k = abs(fft(PSI'))/Nx;                  % Absolute normalized fft
data = log(PSI_k(1:num, 1:J:end));          % Our data

% Prepare the analytical fit
r = sqrt(3)/2;                              % Parameter
c = 37.883;                                 % Parameter
d = 0.5333;                                 % Parameter
e = 1.86;                                   % Parameter

[tg, kg] = meshgrid(t, 1:num-1);            % mesh t and k for 2D calculations
A_k = exp(-kg*r.*sqrt(d+e*(tg-c/r).^2));    % A_k for all k and t
A_0 = sqrt(1 - 2*sum(A_k.^2));              % A_0


% Plot
figure
plot(t, log(A_0), '-')
hold all
plot(t, log(A_k), '-')
hold all
plot(t(1:J:end), data, '+'); grid on;
ylim([-50, 5]);
    
end