function [error, E, ke, pe] = energy(psi, t, k2, Nx, gamma, dt)
% FUNCTION finds the change in energy dE throughout the whole simulation
% INPUT:
%       psi: a matrix with the wave function evolution in space and time.
%            rows are time, columns are space.
%       t: time array
%       k2: square of wave number
%       Nx: number of Fourier modes/spatial nodes
%       gamma: strength of the nonlinearity
%       dt: temporal separation
% OUTPUT:
%       error: integrated energy error
%       E: array with energy at every time step
%       ke: array with kinetic energy at every time step
%       pe: array with potential energy at every time step


pe = zeros(1, length(t));                           % Preallocate pe
ke = zeros(1, length(t));                           % Preallocate ke

for i = 1:length(t)                                 % loop over time
    psi2 = abs(psi(i, :)').*abs(psi(i, :)');        % square psi
    pe(i) = sum(0.5*gamma*psi2.*psi2)/sum(psi2);    % calculate PE

    psi_k = fftshift(fft(psi(i, :)'))/Nx;           % Calculate normalized FFT
    psi_k2 = abs(psi_k).*abs(psi_k);                % Square normalized FFT
    ke(i) =  sum(0.5*k2.*psi_k2)/sum(psi_k2);       % Calculate ke
end

E = ke + pe;                                        % Sum PE and KE to get E

E0 = E(1);                                          % Initial energy
dE = E-E0;                                          % Energy error

figure
plot(t, dE, '-b', 'LineWidth', 2)                   % Plot
xlabel('t'); ylabel('E-E0');                        % Legend
title(sprintf('dt = %0.4f', dt))                    % Title

error = sum(dE)*dt;                                 % integrated energy error
end

