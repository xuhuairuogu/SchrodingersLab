function psi = T2(psi, dt, k2, gamma)
% phi2
%   This function calculates one time step using a second order split step
%   algorithm. 
%   OUTPUT:
%       psi: evolution after one time step.
%   INPUT:
%       psi: function in spatial domain (not fft'd).
%       dt:  temporal separation between nodes.
%       k2:  square of wave numbers from -Nx/2 to Nx/2+1 * whatever. Do not
%            fft shift this.
%       gamma: nonlinearity.
%
%   FORTRAN CODE:
%       epot=agp*cdabs(psi)**2
%       phi = cdexp( -im * 0.5d0*dt * epot )*psi
%       call fft (phi, m, 0 )
%       phi =cdexp( -im * dt * k2/2 )*phi
%       call fft (phi, m, 1 )
%       epot=agp*cdabs(phi)**2
%       psi = cdexp( -im * 0.5d0*dt * epot )*phi

psi2 = gamma*abs(psi).^2;                % Calculate potential
psi = exp(-1i * dt/2 * psi2).*psi;       % Nonlinear calculation
psi = fftshift(fft(psi));                % FFT
psi = exp(-1i * dt * k2/2).*psi;         % Linear calculation
psi = ifft(fftshift(psi));               % Inverse FFT
psi2 = gamma*abs(psi).^2;                % Calculate potential
psi = exp(-1i * dt/2 * psi2).*psi;       % nonlinear calculation

end

