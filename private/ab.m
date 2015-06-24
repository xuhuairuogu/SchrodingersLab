function [wave_2] = ab(full, x, t, row_index, column_index, peak)
% Prepare for spatial profile
wave_1 = full(:,column_index);
x_shift = x(row_index);
t_shift = t(column_index);
%t = linspace(-pi, pi, 1000);
t_1 = 0;
x_1 = x;

a = 3/8;
b = sqrt(8*a*(1-2*a));
omega = 2*sqrt(1-2*a);
psi_1 = (1 + (2*(1-2*a)*cosh(b*t_1) + 1i*b*sinh(b*t_1))./(sqrt(2*a)*cos(omega*(x_1-x_shift))-cosh(b*t_1))).*exp(1i*t_1);

figure
plot(x_1, abs(psi_1).^2, '-r', 'LineWidth', 2); grid on;
hold on
plot(x_1, abs(wave_1).^2, 'b+', 'MarkerFaceColor', 'b');
xlim([-pi, pi]);
tpeak = max(abs(psi_1).^2);
str = sprintf('x=%.3f, t=%.3f, peak=%.3f, th-peak=%.3f.', x_shift, t_shift, peak, tpeak);
title(str);

% prepare for temporal profile
wave_2 = full(row_index,:);
wave_2 = wave_2(1:8:end);
t = t(1:8:end);
t_2 = t-t_shift;
x_2 = 0;
psi_2 = (1 + (2*(1-2*a)*cosh(b*t_2) + 1i*b*sinh(b*t_2))./(sqrt(2*a)*cos(omega*x_2)-cosh(b*t_2))).*exp(1i*t_2);
figure
plot(t, abs(psi_2).^2, '-r', 'LineWidth', 2); grid on;
hold on
plot(t, abs(wave_2).^2, 'b+', 'MarkerFaceColor', 'b');
xlim([35 50]);
tpeak = max(abs(psi_2).^2);
str = sprintf('x=%.3f, t=%.3f, peak=%.3f, th-peak=%.3f.', x_shift, t_shift, peak, tpeak);
title(str);
xlabel('t');
ylabel('|\psi|^2');