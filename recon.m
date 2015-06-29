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

function [psi] = recon(c, d, e, Nx, Nt, tMax)
% RECON: 
%           Reproduces the rogue wave from a specific simulation using some
%           paramters c, d, e.
% INPUT:
%           c: parameter 1, for controlling temporal position of peak
%           d: parameter 2, for controlling height of peak
%           e: parameter 3, for controlling width of peak
%           Nx: number of fourier modes
%           tMax: maximum time to run simulation. 50 is advised. Model
%                 doesn't work after a certain limit past the rogue wave
%                 peak
%           Nt: number of elements of time array
% OUTPUT:
%           psi: full reproduced spatiotemporal evolution
%       

% Calculations to determine minimum value of A_0
t = linspace(0, tMax, Nt);                % Allocate time
r = sqrt(3)/2;                            % Slope
%c = 37.883;                              % Intercept
%d = 0.5333;                              % Parameter
%e = 1.86;
%Nx = 128;                                % Number of fourier modes
k = 1:Nx/2;                               % Wave number
[tg, kg] = meshgrid(t, k);                % mesh t and k for 2D calculations
A_k = exp(-kg*r.*sqrt(d+e*(tg-c/r).^2));  % A_k
A_0 = sqrt(1 - 2*sum(A_k.^2));            % A_0

[~, ~, V] = find(real(A_0));              % Find where min A_0 occurs
index = (real(A_0) == min(V));            % Index of minimum value                  

% Recreate the whole wave
x = linspace(-pi, pi, 101);
psi = zeros(length(t), length(x));        % Final result will be saved here
A_k = [A_0; A_k];
for m = 1:length(t)                       % Outer loop over time
    for xl = x                            % Second loop over space
        total = 0;                        % psi(x = xl, t = m) accumulated here 
        for k = (0:Nx/2)                  % Loop over k to compute psi(xl, m)
            if k == 0                     % k = 0 gets special treatment
                f = 1;
            else
                f = 2;                    % For all values of k
            end                             
            total = total + ...           % psi(x=xl, t=m) = A_0 + ...           
                    f*A_k(k+1, m)...      % 2*sum_k(A_k(m)*cos(k*xl))
                    *cos(k*xl);
        end
        i = (x == xl);                    % Find index of current
        psi(m, i) = total;                % Save psi(xl, m).
    end
    disp([m length(t)]);                  % Show progress                         
end

% Plot surface
figure
surf(x,t,abs(psi).^2, 'EdgeColor', 'none')      % Plot
hold on
ylim([0, 50])                                   % Limit
xlim([x(1), x(length(x))])                      % Limit
colorbar('eastoutside')                         % Colorbar
ylabel('t'); xlabel('x'); zlabel('|\psi|^2')    % Labels  

ab(psi, x, t, max(max(abs(psi).^2)), 3/8);      % Compare to 

end