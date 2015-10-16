% Copyright 2015, Omar Ashour.
% This sourcecode is available from <https://github.com/oashour/HighNLSE/>
%
% This file is part of HighNLSE.
% 
% HighNLSE is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% HighNLSE is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with HighNLSE.  If not, see <http://www.gnu.org/licenses/>.

function psi = T4_NS(psi, dt, k2, gamma)

s = 2^(1/3);
os = 1/(2-s);

ft = os;
bt = -s*os;

psi = T2(psi, ft*dt, k2, gamma);
psi = T2(psi, bt*dt, k2, gamma);
psi = T2(psi, ft*dt, k2, gamma);

end
