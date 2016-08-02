function psi = T6S(psi, dt, k2, V, x, mult)
% T6S:
%   This function calculates one time step using a sixth order split step
%   symplectic algorithm. 
%   OUTPUT:
%       psi: evolution after one time step.
%   INPUT:
%       psi: function in spatial domain (not fft'd).
%       dt:  temporal separation between nodes.
%       k2:  square of wave numbers. Do not fft shift this.
%       V:   potential
s = 2^(1/5);
os = 1/(2-s);

ft = os;                % Forward Step
bt = -s*os;             % Backward Step

psi = T4S(psi, ft*dt, k2, V, x, mult);
psi = T4S(psi, bt*dt, k2, V, x, mult);
psi = T4S(psi, ft*dt, k2, V, x, mult);