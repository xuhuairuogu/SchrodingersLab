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

function [A_0, t] = recon(Nx, Nt, tMax, peakTime)

close all

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
[tg, kg] = meshgrid(t, 1:Nx/2);

q = 3/8;
a = sqrt(8*q*(1-2*q));
alpha = sqrt(2*q)./cosh(a*(tg-peakTime));
alpha_2 = sqrt(2*q)./cosh(a*(t-peakTime));

A_0 = 1 - ((2*(1-2*q)+1i*a*tanh(a*(t-peakTime)))./sqrt(1-alpha_2.^2));
A_k = -(2*(1-2*q)+1i*a*tanh(a*(tg-peakTime)))./sqrt(1-alpha.^2).*((1-sqrt(1-alpha.^2))./alpha).^abs(kg);

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
%psi = [psi; flipud(psi)];
%t = [t t+tMax];
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

% c =37.88
% d = (log(3)/sqrt(3))^2
% e = 0.19