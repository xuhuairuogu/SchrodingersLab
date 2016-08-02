function [psi] = T2(psi, dt, k2, V, x, mult)
% T2:
%   This function calculates one time step using a second order split step
%   algorithm. 
%   OUTPUT:
%       psi: evolution after one time step.
%   INPUT:
%       psi: function in spatial domain (not fft'd).
%       dt:  temporal separation between nodes.
%       k2:  square of wave numbers. Do not fft shift this.
%       V:   potential
Nx = length(x);

pot = V(psi, x);                        % Calculate potential
psi = exp(-1i * pot .* dt/2).*psi;        % Nonlinear calculation

psi = fftshift(fft(psi));                % FFT
psi = exp(-1i * dt * k2/2).*psi;         % Linear calculation
psi = ifft(fftshift(psi));               % Inverse FFT

pot = V(psi, x);                         % Calculate potential
psi = exp(-1i * pot .* dt/2).*psi;       % Nonlinear calculation

% See documentation for explanation of what mult does.
% You will almost always want to set mult to 1 unless you are investigating
% Nonlinear talbot carpets formed by AB's. For example, let's say your mult
% is set to 3, i.e. box size is 3*periodic length. This will result in the
% non-triplet modes growing and ruining your carpet, so this sets them to 0
% and kills their growth artificially. Pretty much useful for nothing but
% generating ugly dress patterns for your grandmother.
if mult > 1                             
    psi = fft(psi);              
    for i = 2:Nx/2+1;
        if(mod(i-1, mult) ~= 0)
            psi(i) = 0;
            psi(Nx - i  + 2) = 0;
        end
    end
    
    psi = ifft(psi);
end

end