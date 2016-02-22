close all; clear variables;

L = 100;
Nx = 1024;
dx = L/Nx;                             % Spatial step size
xo = (-Nx/2:1:Nx/2-1)'*dx;               % Spatial grid points
for M=1:1

psi_i = -0.43 + 0.44*sech(xo/M);
psi_r = -0.84-0.19*sech(xo/M).^2;
psi_0 = psi_r' + 1i*psi_i';

plot(xo, psi_r, '-r', 'LineWidth', 1.5); hold on; plot(xo, psi_i, '-b', 'LineWidth', 1.5); grid on; legend('Real', 'Imag', 0);
title(sprintf('M = %d', M));
hold off;
drawnow;
end