function [] = ab(wave, dx, dt, peak, x)
%t = linspace(-pi, pi, 1000);
t = 0;
a = 3/8;

b = sqrt(8*a*(1-2*a));
omega = 2*sqrt(1-2*a);
psi = 1 + (2*(1-2*a))./(sqrt(2*a)*cos(omega*(x-dx))-1);

figure
plot(x, abs(psi).^2, '-r', 'LineWidth', 2); grid on;
hold on
plot(x, abs(wave).^2, 'b+', 'MarkerFaceColor', 'b');
xlim([-pi, pi]);
disp(max(abs(psi).^2));
str = sprintf('x: %.3f, t: %.3f, peak: %.3f.', dx, dt, peak);
title(str);