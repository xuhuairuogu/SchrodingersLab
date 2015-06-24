function [E ke pe] = energy(u, k2, Nx, gamma)

%          sum_i [(1/2)g*|psi_i|^2] |psi_i|^2
%E_pot= ----------------------------------------
%           sum_i |psi_i|^2

%           sum_i [(1/2)k_i^2] |psik_i|^2
% E_kin= ----------------------------------------
%            sum_i |psik_i|^2

u2 = abs(u).*abs(u);
pe = sum(0.5*gamma*u2.*u2)/sum(u2);
uk = fftshift(fft(u))/Nx;
uk2 = abs(uk).*abs(uk);
ke =  sum(0.5*k2.*uk2)/sum(uk2);
E = ke + pe;
end

