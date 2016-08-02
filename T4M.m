function psi = T4M(psi, dt, k2, V, x, mult)
% T4M:
%   This function calculates one time step using an fourth order split step
%   multi-product algorithm. 
%   OUTPUT:
%       psi: evolution after one time step.
%   INPUT:
%       psi: function in spatial domain (not fft'd).
%       dt:  temporal separation between nodes.
%       k2:  square of wave numbers. Do not fft shift this.
%       V:   potential (function handle)

psi1 = T2(psi, dt/3, k2, V, x, mult);    
psi1 = T2(psi1, dt/3, k2, V, x, mult);  
psi1 = T2(psi1, dt/3, k2, V, x, mult); 
psi1 = 81/40 * psi1;

psi2 = T2(psi, dt/2, k2, V, x, mult); 
psi2 = T2(psi2, dt/2, k2, V, x, mult); 
psi2 = -16/15*psi2;

psi3 = T2(psi, dt, k2, V, x, mult); 
psi3 = 1/24*psi3;

psi = psi1 + psi2 + psi3;

end