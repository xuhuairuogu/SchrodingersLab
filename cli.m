% MAIN CLI SOLVER FOR SchrodingerMAT
clear variables
close all

Nx = 2^8;                   % Number of modes
dt = 0.001;  % Temporal spacing delta t
Tmax = 100;     % Max time for simulation
Nt = Tmax/dt;
mult = 1;

V = @(psi, x) (-1*abs(psi).^2);          % Potential
method = 'T2';                        % Algorithm to use, this means 4th order algorithm with multi-product integrators

gamma0 = 20; % Other parameters for Absorbing BC
alpha = 1;   % Other parameter for absorbing BC
gamma = @(x, gamma_0, alpha, L) gamma_0*(sech(alpha*(x-L/2)).^2 + sech(alpha*(x+L/2)).^2);        % Absorbing boundary conditions function
absorption = struct('useAbsorbingBC', 0, 'gamma', gamma, 'gamma0', gamma0, 'alpha', alpha);       %absorbing BC structure

f1 = figure;
h1 = axes();
f2 = figure;
h2 = axes();
f2.Visible = 'off';
f3 = figure;
h3 = axes();
f3.Visible = 'off';
handles = struct('axes1', h1, 'axes2', h2, 'axes3', h3);

f2.Position = [100, 100, 560*2, 420*2]; % VIDEO
f3.Position = [100, 100, 560*2, 420*2]; % VIDEO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USUAL INITIAL COND +++ DT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2 --> 0.42, 5
%3 --> 0.46, 5
%4 --> 0.47395, 15
%a = 0.5*(1-1/(order+1)^2);

% order = 2;
% a = 0.42;
% L = pi/sqrt(1-2*a);
% Omega = 2*sqrt(1-2*a);
% dx = L/Nx;                             % Spatial step size
% x = (-Nx/2:1:Nx/2-1)'*dx;
% approximate = 1;

%[psi_dt, x_dt] = calcDarboux(order, a, Nx, 0, Tmax/2, 0);
%psi_0 = psi_dt';
%fitP = genCoeff(x_dt, psi_dt, a, order, 0, 0);
%psi_r = fitP{1}(x_dt);
%psi_i = fitP{2}(x_dt);
%coeff = coeffvalues(fitP{1}) + 1i*coeffvalues(fitP{2});
%disp(coeff);

%psi_0 = psi_r + 1i*psi_i;
%psi_0 = psi_dt';
% if approximate
%     A1 = fix(abs(coeff(2))*1e2)/1e2;
%     A2 = fix(abs(coeff(3))*1e2)/1e2;
% else
%     A1 = abs(coeff(2)); A2 = abs(coeff(3));
% end
% A0 = sqrt(1-2*A1^2-2*A2^2);
% phi_0 = angle(coeff(1));
% phi_1 = angle(coeff(2));
% phi_2 = angle(coeff(3));
% A1 = A1*exp(1i*(phi_1 - phi_0));
% A2 = A2*exp(1i*(phi_2 - phi_0));
% psi_0 = A0 + 2*A1*cos(Omega*x) + 2*A2*cos(2*Omega*x);
% norm = dx*(psi_0')*(psi_0)/L;
% disp(['Norm is: ', num2str(norm)]);
% disp([A0, A1, A2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Darboux Transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%psi_0 = calcDarboux(nu, Nx, Nt, 6.3*2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Triplet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L = 20;
% beta = 0;
% gamma = 0;
% t0 = 0;
% psi_0 = triplet(beta, gamma, Nx, 200, L, t0, t0+2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARTIFICIALLY SMOOTH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L = 100;
% M = 1;
% dx = L/Nx;                             % Spatial step size
% xo = (-Nx/2:1:Nx/2-1)'*dx;               % Spatial grid points
% psi_i = -0.43 + 0.44*sech(xo/M);
% psi_r = -0.84-0.19*sech(xo/M).^2;
% psi_0 = psi_r' + 1i*psi_i';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARTIFICAL 2 AB's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%a1 = 0.375;
% nu = 2;
% a1= 0.5*(1-1/nu^2);
% L = 5*pi/sqrt(1-2*a1);
% v1 = 0;
% psi_0 = analytical_ab(v1, a1, Nx, 100, L, 0, 5);
% psi_02 = analytical_ab(-v1, a1, Nx/2, 100, L/2, 0, 5);
% psi_0 = [psi_01, psi_02];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_0 = 0.375;
da = 0.005;
a_f = 0.49;

v1 = VideoWriter('results/intensity.avi');
v1.FrameRate = 1;
v2 = VideoWriter('results/spectrum.avi');
v2.FrameRate = 1;
open(v1);
open(v2);

tic
for a = a_0:da:a_f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYTICAL AB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Omega = 2*sqrt(1-2*a);
% mult = 1;     % KEEP AT ONE EXCEPT WITH NL TALBOT CARPETS IN WHICH CASE UNCOMMENT STUFF FROM T2
% L = mult*pi/sqrt(1-2*a);
% dx = L/Nx;     
% x = (-Nx/2:1:Nx/2-1)'*dx;  
% t_s = 0;
% x_s = x;
% lambda = Omega*sqrt(1-Omega^2/4);
% psi_s = ((1-4*a)*cosh(lambda*t_s)+sqrt(2*a)*cos(Omega*x_s)+1i*lambda*sinh(lambda*t_s))./(sqrt(2*a)*cos(Omega*x_s)-cosh(lambda*t_s));
% psi_0 = psi_s;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USUAL INITIAL COND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A1 = 1e-2;
A2 = 0;
A0 = sqrt(1-2*A1^2-2*A2^2);

Omega = 2*sqrt(1-2*a);
L = pi/sqrt(1-2*a);
dx = L/Nx;                             % Spatial step size
x = (-Nx/2:1:Nx/2-1)'*dx;   
psi_0 = A0 + 2*A1*cos(Omega*x) + 2*A2*cos(2*Omega*x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[PSI, x, t] = solve(dt, Nx, Tmax, L, mult, V, psi_0, absorption, method);

densityPlot(abs(PSI).^2, x, t, ceil(Nt/1000), ceil(Nx/256), h2); colormap('jet'); caxis([0, 10]);
PSI_k = log(abs(fft(PSI'))/length(PSI(1, :)));
fourierPlot(PSI_k.', t, 7, 1, h3); legend off; ylim([-25, 5]);
title(h2, sprintf('a = %.4f', a));
title(h3, sprintf('a = %.4f', a));
hold off;


frame1 = getframe(f2);
writeVideo(v1,frame1);
frame2 = getframe(f3);
writeVideo(v2,frame2);
   
end
close(v1);
close(v2); 
toc;

close(f1);