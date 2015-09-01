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

function [spatial, temporal, t_shift] = ab(psi, x, t, peak, a)
% FUNCTION: [spatial temporal] = ab(psi, x, t, peak)
%           Compares a full time evolution of a wave function psi to an Akhmediev
%           breather using supplied value of a
% INPUT:    psi:  Matrix representing wafe function with rows as time and
%                 columns as space.
%           x:    spatial variable used to create psi
%           t:    temporal variable used to create psi
%           peak: actual peak to compare to theoretical peak
%           a:    parameter used for akhmediev breather
% OUTPUT:   spatial: spatial profile array
%           temporal: temporal profile array

% Find row and column of the rogue wave peak
[~, maxIndex] = max(psi(:));                    %
[row, col] = ind2sub(size(psi), maxIndex);      %

x_shift = x(col);                               % x-value of the peak
t_shift = t(row);                               % t-value of the peak

% Spatial profile
spatial = psi(row,:);                           % spatial profile corresponds to 
                                                % a row of psi
t_s = 0;                                        % Use t = 0 for analytical AB
x_s = x - x_shift;                              % Shift space for analytical AB

b = sqrt(8*a*(1-2*a));                          % Paramter
omega = 2*sqrt(1-2*a);                          % Parameter
psi_s = (1 + (2*(1-2*a)*cosh(b*t_s) + 1i*b*sinh(b*t_s))./ ... % Analytical
        (sqrt(2*a)*cos(omega*x_s)-cosh(b*t_s))).*exp(1i*t_s); % Spatial profile

% Temporal profile
temporal = psi(:,col);                          % temporal profile corresponds 
                                                % to column of psi
t_t = t - t_shift;                              % Shift time
x_t = 0;                                        % Use x = 0 for analytical AB
psi_t = (1 + (2*(1-2*a)*cosh(b*t_t) + 1i*b*sinh(b*t_t))./ ... % Analytical
        (sqrt(2*a)*cos(omega*x_t)-cosh(b*t_t))).*exp(1i*t_t); % Temporal profile

% Plot the results
figure
plot(x, abs(psi_s).^2, '-r', 'LineWidth', 2); % Analytical spatial
grid on; hold on;                             % Grid
plot(x, abs(spatial).^2, 'b+', ...            % Actual spatial
     'MarkerFaceColor', 'b', 'MarkerSize', 8);
xlim([-pi, pi]);                              % Limits
xlabel('x'); ylabel('|\psi|^2|');             % Labels
%title(sprintf('x=%.3f, t=%.3f, peak=%.3f, th-peak=%.3f.', ... % Title
%              x_shift, t_shift, peak, max(abs(psi_s).^2)));   %
legend('Analytical', 'Calculation', 0); 
          
figure
plot(t, abs(psi_t).^2, '-r', 'LineWidth', 2); % Analytical temporal
grid on; hold on;                             % Grid
plot(t(1:10:end), abs(temporal(1:10:end)).^2, 'b+', ...           % Actual temporal
     'MarkerFaceColor', 'b', 'MarkerSize', 8);
xlim([t_shift-5, t_shift+5]);                 % Limits
xlabel('t'); ylabel('|\psi|^2|');             % Labels
%title(sprintf('x=%.3f, t=%.3f, peak=%.3f, th-peak=%.3f.', ... % Title
%              x_shift, t_shift, peak, max(abs(psi_s).^2)));   %
legend('Analytical', 'Calculation', 0);          
          
disp(abs(spatial(1)).^2);
end