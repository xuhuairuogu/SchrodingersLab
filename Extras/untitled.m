close all;
clear variables;

Nx = 10240;
Nz = 500;
Lx = 2*4.412;
Lz = 5;
dx = Lx/Nx;                             % Spatial step size
xo = (-Nx/2:1:Nx/2-1)'*dx;           
dz = Lz/Nz;                             % Spatial step size
zo = (-Nz/2:1:Nz/2-1)'*dz;
[x, z] = meshgrid(xo, zo);

k = 1.2;
psi = AB_je(xo, 0, k);

integrand = @(m, theta) 1./sqrt(1-m^2*sin(theta).^2);
if k > 1
    m = sqrt((k-1)/(2*k));
    L = 4/sqrt(2*k)*integral(@(theta) integrand(m, theta), 0, pi/2,'AbsTol',1e-16,'RelTol',1e-10);
elseif k <= 1
    m = sqrt((1-k)/(1+k));
    L = 4/sqrt(1+k)*integral(@(theta) integrand(m, theta), 0, pi/2,'AbsTol',1e-16,'RelTol',1e-10);
end
disp(L);
%  [~, I] = max(abs(psi(length(z)/2+1, :)).^2);
%  disp(0-xo(I));

% surf(xo, zo, abs(psi).^2, 'EdgeColor', 'none');
% view([0 0 90]);
% xlabel('x');
% ylabel('t');
% colorbar;
% colormap('jet');

% figure
% plot(xo, abs(psi(length(zo)/2-1, :)).^2);

figure
plot(xo, abs(psi).^2);