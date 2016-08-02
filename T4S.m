function psi = T4S(psi, dt, k2, V, x, mult)
% T4S:
%   This function calculates one time step using a fourth order split step
%   symplectic algorithm. 
%   OUTPUT:
%       psi: evolution after one time step.
%   INPUT:
%       psi: function in spatial domain (not fft'd).
%       dt:  temporal separation between nodes.
%       k2:  square of wave numbers. Do not fft shift this.
%       V:   potential
s = 2^(1/3);
os = 1/(2-s);

ft = os;                % Forward factor
bt = -s*os;             % Backward factor
 
psi = T2(psi, ft*dt, k2, V, x, mult);
psi = T2(psi, bt*dt, k2, V, x, mult);
psi = T2(psi, ft*dt, k2, V, x, mult);
end
