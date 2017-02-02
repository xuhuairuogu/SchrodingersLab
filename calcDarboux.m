function [PSI, xo, to, L] = calcDarboux(n, a, R, tmax, g, seed, xs, ts, Nx, Nt, mode, Lk, mult)
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
% This needs to be rewritten to use "l" instead of a.
if length(a) == 1 && R == 2 && strcmp(seed, 'Breather') && n >= 2;        % Maximal intensity family mode
    for k = 1:n
        a(k) = k^2*(a(1)-1/2)+1/2;
    end
elseif length(a) == 1 && R == 2 && strcmp(seed, 'Dn') && n >= 2;        % Maximal intensity cnoidal family mode
    v = imag(a);
    for k = 1:n
        P = (g.^4.*k.^2+8.*((-2)+g.^2).*(-1+k.^2).*v.^2+16.*k.^2.*v.^4).^(1/2);
        nu(k) = (1/4).*2.^(-1/2).*v.^(-1).*(P.^2+(P.^4+(-64).*g.^4.*v.^4).^(1/2)).^(1/2);
    end
    a = 1i*nu;
elseif length(a) == 1 && R ~= 2;    % Second order RW where Omega2=R*Omega1
    if n > 2
        uiwait(warndlg('Ratio not equal to 2 and order > 2. Feature not yet implemented')); 
        return
    end 
    a2 = R^2*(a(1) - 1/2) + 1/2;
    a = [a(1), a2];
end

% Determine whether in eigenvalue or 'a' mode.
if ~imag(a) % If 'a' mode, calculate purely imag eigenvalues
    l = 1i*sqrt(2*a); 
else % Otherwise, eigenvalues already entered
    l = a;
end

if strcmp(seed, 'Breather') || strcmp(seed, 'Soliton');
    kappa = 2*sqrt(1+l.^2);      % Principal wave number
    chi = 0.5*acos(kappa/2);
elseif strcmp(seed, 'Dn')
    kappa = sqrt(1+(l - g^2/4./l).^2);  % Half the principal wave number
    chi = 0.5*acos(kappa);
elseif strcmp(seed, 'Cn')
    kappa = g*sqrt(1+1/g^2*(l - 1/4./l).^2); 
    chi  = 0.5*acos(kappa/g); 
end

% Denominator of ratio between Omega2 and Omega1 determines box size. This
% has not been published yet.
[~, D] = rat(R); 

% Decide box size calculation mode
switch mode
    case 'manual' % Manual mode
        L = Lk;
    case 'periodic' % Periodic mode
        switch lower(seed)
            case 'breather'
                L = 2*D*mult*pi/kappa(1);     % Periodic length
            case 'dn'
                L = D*mult*pi/kappa(1);
            otherwise
                error('Mode set to periodic for a non-periodic solution');
        end
end

% Determine box type. Essentially, if Nt > 0, then we are calculating the
% usual t array. Same applied for x. If either of those are zero, then we
% are interested in the value of the wave only at t=tMax or x = xo. i.e.
% either spatial or temporal profile of the breather/RW at a specific
% point. If both are 0, we are calculting only 1 point using the DT. This
% is useful for comparing the DT peak to the analytical formula we
% published in PLA.
if Nt > 0 
    dt = tmax/Nt;
    to = (-Nt/2+1:1:Nt/2)'*2*dt;                                % Propagation variable
else 
    to = tmax;
end

if Nx > 0
    dx = L/Nx;                                                  % Spatial step size
    xo = (-Nx/2+1:1:Nx/2)'*dx; 
else
    xo = 0;
end
% Begin recursion
psi = cell(n, 1);
[~, ~, psi] = calc_rs(n, 1, xo, to, g, seed, xs, ts, l, kappa, chi, psi);       

% Output the actual solution we are interested in
PSI = psi{n};    

function [rf, sf, psi] = calc_rs(n, p, x, t, g, seed, xs, ts, l, kappa, chi, psi)
    % calc_rs: Recursive function for calculation of the Darboux
    %          Transformation for investigation analytical solutions of the
    %          nonlinear Schrodinger equation.
    % INPUT: n: order of the current step. 
    %        p: enumeration counter.
    %        x: transverse variable
    %        t: evolution variable
    %        xs: array of x-shifts
    %        ts: array of t-shifts
    %         l: complex eigenvalue array
    %        kappa: array of wavenumbers
    %        psi: cell array that holds each step of the transform
    % For more details, see: CITE AKHMEDIEV ET AL.
    
    if n == 1 % Base case
        
        % Complex lambda stokes' seeding 
        if strcmp(seed, 'Breather')
            [x,t] = meshgrid(x,t);
            dd = l(p)*kappa(p);
            A = +chi(p) + 0.5*(kappa(p)*(x-xs(p))+dd*(t-ts(p))) - pi/4;
            B = -chi(p) + 0.5*(kappa(p)*(x-xs(p))+dd*(t-ts(p))) - pi/4;
            rf = 2*1i*exp(-1i*t/2).*sin(A);
            sf = 2*exp(1i*t/2).*cos(B);
        elseif strcmp(seed, 'Soliton')
            [x,t] = meshgrid(x,t);
            rf = exp(1i*l(p)*(x-xs(p)) + 1i*l(p)^2*(t-ts(p)) - 1i*pi/4);
            sf = exp(-1i*l(p)*(x-xs(p)) - 1i*l(p)^2*(t-ts(p)) + 1i*pi/4);
        elseif strcmp(seed, 'Dn')
            [~,~,DN] = ellipj(x,g^2);
            A = exp(-1i*pi/4); B = exp(+1i*pi/4);
            Nx = length(x); dx = x(2) - x(1); iz = find(x == 0);
            a = zeros(length(x), length(t)); b = zeros(length(x), length(t));
            a(iz, :) =  A*exp(1i*(chi(p) + kappa(p)*l(p)*(t-ts(p)))) ...
                -B*exp(-1i*(chi(p) + kappa(p)*l(p)*(t-ts(p))));
            b(iz, :) =  A*exp(1i*(-chi(p) + kappa(p)*l(p)*(t-ts(p)))) ...
                +B*exp(-1i*(-chi(p) + kappa(p)*l(p)*(t-ts(p))));

            for i=1:Nx/2
                a(iz+i, :) = a(iz+i-1,:) + dx*(1i*l(p)*a(iz+i-1,:) + ...
                          1i*b(iz+i-1,:)*DN(iz+i-1));
                b(iz+i, :) = b(iz+i-1,:) + dx*(-1i*l(p)*b(iz+i-1,:) + ...
                          1i*a(iz+i-1,:)*DN(iz+i-1));
                if i~=Nx/2
                    a(iz-i,:) = a(iz+i,:);
                    b(iz-i,:) = b(iz+i,:);
                end
            end
            [x,t] = meshgrid(x,t);
            %x = x'; t = t';
            rf = a.'.*exp(+1i*t/4*(g^2-2));
            sf = b.'.*exp(-1i*t/4*(g^2-2));
        elseif strcmp(seed, 'Cn')
            [~,CN,~] = ellipj(x,g^2);
            A = exp(-1i*pi/4); B = exp(+1i*pi/4);
            Nx = length(x); dx = x(2) - x(1); iz = find(x == 0);
            a = zeros(length(x), length(t)); b = zeros(length(x), length(t));
            a(iz, :) =  A*exp(1i*(chi(p) + kappa(p)*l(p)*(t-ts(p)))) ...
                -B*exp(-1i*(chi(p) + kappa(p)*l(p)*(t-ts(p))));
            b(iz, :) =  A*exp(1i*(-chi(p) + kappa(p)*l(p)*(t-ts(p)))) ...
                +B*exp(-1i*(-chi(p) + kappa(p)*l(p)*(t-ts(p))));

            for i=1:Nx/2
                a(iz+i, :) = a(iz+i-1,:) + dx*(1i*l(p)*a(iz+i-1,:) + ...
                          1i*b(iz+i-1,:)*g*CN(iz+i-1));
                b(iz+i, :) = b(iz+i-1,:) + dx*(-1i*l(p)*b(iz+i-1,:) + ...
                          1i*a(iz+i-1,:)*g*CN(iz+i-1));
                if i~=Nx/2
                    a(iz-i,:) = a(iz+i,:);
                    b(iz-i,:) = b(iz+i,:);
                end
            end
            [x,t] = meshgrid(x,t);
            rf = a.'.*exp(-1i*t/4*(2*g^2-1));
            sf = b.'.*exp(+1i*t/4*(2*g^2-1));
        else
            error('Unknown Seed.');
        end
 
        if (p == 1)
            if strcmp(seed, 'Breather')
                psi_0 = exp(1i*t);
            elseif strcmp(seed, 'Soliton')
                psi_0 = 0;
            elseif strcmp(seed, 'Dn')
                [~,~,DN] = ellipj(x,g^2);
                psi_0 = DN.*exp(1i*(t-ts(p))*(1-g^2/2));
            elseif strcmp(seed, 'Cn');
                [~,CN,~] = ellipj(x,g^2);
                psi_0 = g*CN.*exp(1i*(t-ts(p))*(g^2-1/2));
            else
                error('Unknown Seed');
            end
            psi{n} = psi_0 + (2*(conj(l(n)) - l(n))*sf.*conj(rf))./(abs(rf).^2 + abs(sf).^2);
        end
    else % Recursion
        [r1, s1, psi] = calc_rs(n-1, 1,   x, t, g, seed, xs, ts, l, kappa, chi, psi);
        [r2, s2, psi] = calc_rs(n-1, p+1, x, t, g, seed, xs, ts, l, kappa, chi, psi);

        rf = ((conj(l(n-1) ) -      l(n-1) ) * conj(s1).*r1.*s2  ...
             +(     l(p+n-1) -      l(n-1) ) * abs(r1).^2.*r2    ...
             +(     l(p+n-1) - conj(l(n-1))) * abs(s1).^2.*r2)   ...
             ./(abs(r1).^2 + abs(s1).^2);
        sf = ((conj(l(n-1) ) -      l(n-1) ) * s1.*conj(r1).*r2  ...
             +(     l(p+n-1) -      l(n-1) ) * abs(s1).^2.*s2    ...
             +(     l(p+n-1) - conj(l(n-1))) * abs(r1).^2.*s2)   ...
             ./(abs(r1).^2 + abs(s1).^2);
       if (p == 1)
            psi{n} = psi{n-1} + (2*(conj(l(n)) - l(n))*sf.*conj(rf))./(abs(rf).^2 + abs(sf).^2);
        end
    end