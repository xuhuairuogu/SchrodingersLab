function [psi] = T1(psi, dt, k2, V, x)
% T1:
%   This function calculates one time step using a first order split step
%   algorithm. 
%   OUTPUT:
%       psi: evolution after one time step.
%   INPUT:
%       psi: function in spatial domain (not fft'd).
%       dt:  temporal separation between nodes.
%       k2:  square of wave numbers. Do not fft shift this.
%       V:   potential

pot = V(psi, x);                         % Calculate potential
psi = exp(-1i * pot .* dt).*psi;         % Nonlinear calculation %%% CHANGE 2/2 to dt/2

psi = fftshift(fft(psi));                % FFT
psi = exp(-1i * dt * k2/2).*psi;         % Linear calculation %%%% CHANGE dt*2 to dt
psi = ifft(fftshift(psi));               % Inverse FFT

end