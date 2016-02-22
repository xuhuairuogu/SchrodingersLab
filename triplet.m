function [psi_0, xo] = triplet(beta, gamma, Nx, Nt, L, tmin, tmax)

clear vairbales
close all;
% Their notation
% G2 = 12*(3 - 16*t.^4 - 24*t.^2*(4*x.^2 + 1) - 4*beta*t - 80*x.^4 - 72*x.^2 + 4*gamma*x)
% K2 = 24*(x*(15-16*t.^4 + 24*t.^2-4*beta*t)-8*(4*.t^2+1)*x.^3 - 16*x.^5 + gamma*(2*x.^2-2*t.^2 - 1/2));
% D2 = 64*t.^6 + 48*t.^4*(4*x.^2+1) + 12*t.^2*(3-4*x.^2).^2 + 65*x.^6 + 432*x.^4 + 396*x.^2 + 9 + beta*(beta+4*t*(12*x.^2-4*t.^2+3))+gamma(gamma + 4*x*(12*t.^2 - 4*x.^2 - 9));
dx = L/Nx;                             % Spatial step size
xo = (-Nx/2:1:Nx/2-1)'*dx;   
to = linspace(tmin, tmax, Nt);
[x, t] = meshgrid(xo, to);
G2 = 12*(3 - 16*x.^4 - 24*x.^2.*(4*t.^2 + 1) - 4*beta*x - 80*t.^4 - 72*t.^2 + 4*gamma*t);
K2 = 24*(t.*(15-16*x.^4 + 24*x.^2-4*beta*x)-8*(4*x.^2+1).*t.^3 - 16*t.^5 + gamma*(2*t.^2-2*x.^2 - 1/2));
D2 = 64*x.^6 + 48*x.^4.*(4*t.^2+1) + 12*x.^2.*(3-4*t.^2).^2 + 64*t.^6 + 432*t.^4 + 396*t.^2 + 9 + beta*(beta+4*x.*(12*t.^2-4*x.^2+3))+gamma*(gamma + 4*t.*(12*x.^2 - 4*t.^2 - 9));

psi_2 = (1 + (G2 + 1i*K2)./D2).*exp(1i*t);

figure;
h1 = axes();
densityPlot(abs(psi_2).^2, xo, to, ceil(Nx/512), ceil(Nt/500), h1); colormap('jet'); % Plot result

psi_0 = psi_2(1, :);
% plot(xo, real(psi_0), '-b'); title(sprintf('real, t=%d', t(1))); grid on;%ylim([-1.04, -0.84]); grid on;
% figure
% plot(xo, imag(psi_0), '-r'); title(sprintf('imag, t=%d', t(1))); grid on;%ylim([-0.6, 0.1]); grid on;
% figure;
% plot(xo, abs(psi_0).^2, '-k'); title(sprintf('intensity, t=%d', t(1)));  ylim([0.99, 1.08]); grid on;

maxima = regions(psi_2, xo, to);
% disp(sqrt(3)/4*(abs(gamma)^(1/3)));
% PSI_k = log(abs(fft(psi_2'))/Nx);
% figure;
% h2 = axes();
% fourierPlot(PSI_k', to, 6, h2)

end