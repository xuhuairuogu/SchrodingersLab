function peak = peakPredict(a, order, R, seed, g)

switch lower(seed)
    case 'soliton'
        seed_0 = 0;
    case 'breather',
        seed_0 = 1;
    case 'cn'
        seed_0 = g;
    case 'dn'
        seed_0 = 1;
end

% Determine which mode DT is running in
if length(a) == 1 && R == 2 && strcmp(seed, 'Breather') && order >= 2;        % Maximal intensity family mode
    for k = 1:order
        a(k) = k^2*(a(1)-1/2)+1/2;
    end
elseif length(a) == 1 && R == 2 && strcmp(seed, 'Dn') && order >= 2;        % Maximal intensity family mode
    for k = 1:order
        a(k) = 1./(128.*a(1)).*(32.*a(1) - 16.*a(1).*g.^2 - 32.*a(1).*k.^2 + 64.*a(1).^2.*k.^2 + 16.*a(1).*g.^2.*k.^2 + g.^4.*k.^2 + ...
                 sqrt(-256.*a(1).^2.*g.^4 + (64.*a(1).^2.*k.^2 + g.^4.*k.^2 + 16.*a(1).*(-2 + g.^2).*(-1 + k.^2)).^2));
    end
elseif length(a) == 1 && R ~= 2;    % Second order RW where Omega2=R*Omega1
    if n > 2
        uiwait(warndlg('Ratio not equal to 2 and order > 2. Feature not yet implemented')); 
        return
    end 
    a2 = R^2*(a(1) - 1/2) + 1/2;
    a = [a(1), a2];
end

% Determine whether in eigenvalue of 'a' mode.
if ~imag(a) % If 'a' mode, calculate purely imag eigenvalues
    l = 1i*sqrt(2*a); 
else % Otherwise, eigenvalues already entered
    l = a;
end
peak = seed_0 + 2*sum(imag(l(1:order)));

end
