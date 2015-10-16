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

function [A_0, time] = a_plot(PSI, Nx, t, num)
% FUNCTION: Plots analytical function of r, c, d, e vs actual data
% INPUT:
%       PSI: full spatiotemporal wave function
%       Nx: Number of Fourier modes/spatial nodes
%       t: time array
%       num: number of elements to plot. 5 is recommended.
% See the file recon for more info on the analytical fit.

% Prepare the actual data
J = 200;                                    % Scaling down of number of points
PSI_k = abs(fft(PSI'))/Nx;                  % Absolute normalized fft
data = log(PSI_k(1:num, 1:J:end));          % Our data
time = t;
r = sqrt(3)/2;                               % Parameter
[tg, kg] = meshgrid(time, 1:num-1);            % mesh t and k for 2D calculations
c = 37.08/(sqrt(3)/2);
A_k = exp(r*kg.*(tg-c));
A_0 = sqrt(1 - 2*sum(A_k.^2));              % A_0

fit = log([A_0; A_k]);
% Plot
figure
fittt = plot(time, fit, '-', 'LineWidth', 1.5);
hold on
data = plot(t(1:J:end), data, 'o', 'MarkerSize', 6); grid on;
ylim([-50, 5]);
xlim([0, max(t)]);
xlabel('t');
ylabel('log(|A_k|)');
%legend([data, fittt], 'Data', 'Fit')

% % Extra plot
% k = (1:num-1);
% index_1 = find(t == 10);
% index_2 = find(t == 45);
% time = t(index_1:index_2);
% m = zeros(1, length(time));
% 
% gif = 0;
% if gif
%     figure(1)
%     filename = 'test.gif';
% end
% for i = 1:length(time)
%     index = t == time(i);
%     disp([i length(time)])
%     A_peak = data(2:num, index)';
%  
%     result = polyfit(k(1:2),A_peak(1:2),1);
%     m(i) = result(1);
% 
%     if gif
%         plot(k, A_peak); 
%         xlabel('k'); ylabel('log(|A_k)');
%         xlim([min(k) max(k)]); ylim([-40, 0]);
%         drawnow
%         title(sprintf('t = %f', time(i))); 
% 
%        frame = getframe(1);
%        im = frame2im(frame);
%        [imind,cm] = rgb2ind(im,256);
%        if i == 1;
%            imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0);
%        else
%            imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0);
%        end  
%     end
% end
% figure
% plot(time(1:30:end), m(1:30:end), 'b+');
% xlabel('t'); ylabel('Slope');
% anal = -r*sqrt(d+(time-43.74).^2);
% hold on
% plot(time, anal, 'r-', 'LineWidth', 2);
% legend('Data', 'Analytical', 0);

end