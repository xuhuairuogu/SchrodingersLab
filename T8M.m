function psi = T8M(psi, dt, k2, V, x, mult)
% T8M:
%   This function calculates one time step using an eighth order split step
%   multi-product algorithm. 
%   OUTPUT:
%       psi: evolution after one time step.
%   INPUT:
%       psi: function in spatial domain (not fft'd).
%       dt:  temporal separation between nodes.
%       k2:  square of wave numbers. Do not fft shift this.
%       V:   potential

psi1 = T2(psi, dt/4, k2, V, x, mult);   
psi1 = T2(psi1, dt/4, k2, V, x, mult);
psi1 = T2(psi1, dt/4, k2, V, x, mult);
psi1 = T2(psi1, dt/4, k2, V, x, mult);
psi1 = 1024/315 * psi1;

psi2 = T2(psi, dt/3, k2, V, x, mult); 
psi2 = T2(psi2, dt/3, k2, V, x, mult);
psi2 = T2(psi2, dt/3, k2, V, x, mult);
psi2 = -729/280*psi2;

psi3 = T2(psi, dt/2, k2, V, x, mult); 
psi3 = T2(psi3, dt/2, k2, V, x, mult);
psi3 = 16/45*psi3;

psi4 = T2(psi, dt, k2, V, x, mult);
psi4 = -1/360*psi4;

psi = psi1 + psi2 + psi3 + psi4;

end