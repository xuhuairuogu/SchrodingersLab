function [rf, sf, psi] = calc_rs(n, p, x, t, xs, ts, l, kappa, psi)
    if n == 1
        A =  0.5*(acos(kappa(p)/2) + (x-xs(p))*kappa(p) - pi/2) + ...
          1i*0.5*((t-ts(p))*kappa(p)*sqrt(1-kappa(p)^2/4));
        B =  0.5*(-acos(kappa(p)/2) + (x-xs(p))*kappa(p) - pi/2) + ...
          1i*0.5*((t-ts(p))*kappa(p)*sqrt(1-kappa(p)^2/4));
      
        rf = 2*1i*exp(-1i*t/2).*sin(A);
        sf = 2*exp(1i*t/2).*cos(B);
        if (p == 1)
            psi{n} = exp(1i*t) + (2*(conj(l(n)) - l(n))*sf.*conj(rf))./(abs(rf).^2 + abs(sf).^2);
        end
    else
        [r1, s1, psi] = calc_rs(n-1, 1,   x, t, xs, ts, l, kappa, psi);
        [r2, s2, psi] = calc_rs(n-1, p+1, x, t, xs, ts, l, kappa, psi);

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
end