function psi_0 = analytical_ab(v1, a1, Nx, Nt, L, tmin, tmax)

to = linspace(tmin, tmax, Nt);
xo = linspace(-L/2, L/2, Nx);

[x, t] = meshgrid(xo, to);

Omega1 = 2*sqrt(1-2*a1);
lambda1 = sqrt(8*a1*(1-2*a1));

% Basic parameters
l1i = sqrt(2*a1);

l1 = v1 + 1i*l1i;

kappa1 = 2*sqrt(1+l1^2);

k1r = real(kappa1);
k1i = imag(kappa1);

chi1 = 1/2*acos(kappa1/2);

X1r = real(chi1);
X1i = imag(chi1);

d1 = l1*kappa1;
d1r = real(d1);
d1i = imag(d1);

Ar = X1r + 1/2*(k1r*x + d1r*t) - pi/4;
Ai = X1i + 1/2*(k1i*x+d1i*t);
Br = -X1r + 1/2*(k1r*x + d1r*t) - pi/4;
Bi = -X1i + 1/2*(k1i*x+d1i*t);

D1 = cos(2*Br) - cos(2*Ar) + cosh(2*Ai) + cosh(2*Bi);

psi_1 = (1 + 8*1i*l1i./D1.*cosh(Bi-1i*Br).*sinh(Ai + 1i*Ar)).*exp(1i*t);

%figure;
%h1 = axes();
%densityPlot(abs(psi_1).^2, xo, to, 1, 1, h1); colormap('jet');        % Plot result
psi_0 = psi_1(1, :);
plot(xo, abs(psi_0));

end