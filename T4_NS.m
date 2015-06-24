function psi = T4_NS(psi, dt, k2, gamma)

% schem4(dt)=schem2(ft*dt)schem2(bt*dt)schem2(ft*dt)

% where the "forward time" coefficient ft and the "backward time" coefficient bt are defined by
%       ft=os
%       bt=-s*os
%where
%       os=1.d0/(2.d0-s)
%       s=2.d0**(1.d0/3.d0)

s = 2^(1/3);
os = 1/(2-s);

ft = os;
bt = -s*os;

psi = T2(psi, ft*dt, k2, gamma);
psi = T2(psi, bt*dt, k2, gamma);
psi = T2(psi, ft*dt, k2, gamma);

end
