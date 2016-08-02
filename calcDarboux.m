function [PSI, xo, to] = calcDarboux(n, a, R, tmax, xs, ts, Nx, Nt, mode, Lk, mult)
% calcDarboux: Wrapper function for calculating the Darboux Trasformation
%              for investigating analytical solutions of the nonlinear
%              Schrodinger equation.
% Input:       n: order of the DT
%              a: array of length n, contains either the parameter a or the
%                 complex eigenvalues of the DT.
%              R: Ratio between Omega1 and Omega 2 (see below). Set to
%                 2 to investigate the maximal intensity breather family.
%              tmax: Maximum running time of the calculation
%              xs: array of length n containing the x-shifts
%              ts: array of length n containing the t-shifts
%              Nx: Number of x-nodes
%              Nt: Number of t-nodes
%              mode: box size selection mode, string. Either 'manual' or
%                    'periodic'. See below.
%              Lk: manual box size if manual is selected above.
%              mult: How many periods to display if 'periodic' is selected.

% Determine which mode DT is running in
if length(a) == 1 && R == 2;        % Maximal intensity family mode
    for k = 1:n
        a(k) = k^2*(a(1)-1/2)+1/2;
    end
elseif length(a) == 1 && R ~= 2;    % Second order RW where Omega2=R*Omega1
    if n > 2
        uiwait(warndlg('Ratio not equal to 2 and order > 2. Feature not yet implemented')); 
        return
    end 
    a2 = R^2*(a(1) - 1/2) + 1/2;
    a = [a(1), a2];
end

% Denominator of ratio between Omega2 and Omega1 determines box size. This
% has not been published yet.
[~, D] = rat(R); 

% Decide box size calculation mode
switch mode
    case 'manual' % Manual mode
        L = Lk;
    case 'periodic' % Periodic mode
        if ~imag(a) % If parameter 'a' is used.
            aa = a(1); 
        else % Otherwise, complex eigenvalue form is assumed
            aa = imag(a(1))^2/2; % Calculate parameter 'a' from 'l'
        end  
        L = D*mult*pi/sqrt(1-2*aa);     % Periodic length
end

% Determine whether in eigenvalue of 'a' mode.
if ~imag(a) % If 'a' mode, calculate purely imag eigenvalues
    l = 1i*sqrt(2*a); 
else % Otherwise, eigenvalues already entered
    l = a;
end
kappa = 2*sqrt(1+l.^2);                           % Principal wave number

% Determine box type. Essentially, if Nt > 0, then we are calculating the
% usual t array. Same applied for x. If either of those are zero, then we
% are interested in the value of the wave only at t=tMax or x = xo. i.e.
% either spatial or temporal profile of the breather/RW at a specific
% point. If both are 0, we are calculting only 1 point using the DT. This
% is useful for comparing the DT peak to the analytical formula we
% published in PLA.
if Nt > 0 
    dt = tmax/Nt;
    to = (-Nt/2:1:Nt/2-1)'*2*dt;                                % Propagation variable
else 
    to = tmax;
end

if Nx > 0
    dx = L/Nx;                                                  % Spatial step size
    xo = (-Nx/2:1:Nx/2-1)'*dx; 
else
    xo = 0;
end
[x, t] = meshgrid(xo, to);                                      % Mesh for calculation

% Begin recursion
psi = cell(n, 1);
[~, ~, psi] = calc_rs(n, 1, x, t, xs, ts, l, kappa, psi);       

% Output the actual RW we are interested in
PSI = psi{n};                                        

end