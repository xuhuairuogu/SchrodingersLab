function [A0, A] = genCoeff(x, psi, a, ord, R, plot_res, handle, approx)
%CREATEFITS(X_DT,PSI_R,PSI_I)
%  Create fits for function supplied by DT to be used as initial condition
%  for RW/AB of ord "ord"
%
%  Input:
%       x: x-coordinated used with psi below.
%       psi: array containing complex psi from DT
%       a: a parameter from DT.
%       ord: ord of RW, and thus, ord of Fourier Series fit
%       plot_res: produce result of fit + residuals
%  Output:
%      fitresult : a cell-array of fit objects representing the fits.
%      gof : structure array with goodness-of fit info.

% Initialization.
Omega = 2*sqrt(1-2*a);
psi_i = imag(psi.'); psi_r = real(psi.');

% Initialize arrays to store fits and goodness-of-fit.
fitP = cell(2, 1);
gof = struct( 'sse', cell(2,1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

% Set up fit type and options.
fitString = sprintf('A0 + 2*A1*cos(%.16f*x)', Omega); % First ord fit
for i=2:ord % Higher ord fits
    fitString = strcat(fitString, sprintf(' + 2*A%d*cos(%.16f*%.16f*x)', i, R(i), Omega)); 
end
ft = fittype(fitString, 'independent', 'x', 'dependent', 'y' );
opts = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', 0.5*ones(1,ord+1));
opts.Display = 'Off';

% Fit real part
[xData, yData] = prepareCurveData(x, psi_r);
[fitP{1}, gof(1)] = fit(xData, yData, ft, opts);

% Plot fit with data.
if plot_res
    %figure('Name', 'Fits');
    % plot fit
    subplot(2, 2, 1, 'Parent', handle);
    h = plot(fitP{1}, xData, yData );
    %legend(h, 'Data', 'Fit', 'Orientation', 'Horizontal', 'Location', 'Best');
    legend off;
    xlabel('x'); ylabel('\psi_r'); grid on;
    % plot residuals
    subplot(2, 2, 2, 'Parent', handle);
    h = plot(fitP{1}, xData, yData, 'Residuals');
    xlabel('x'); ylabel('\psi_r'); grid on;
    %legend(h, 'Residuals', 'Zero-line', 'Orientation', 'Horizontal', 'Location', 'Best');
    legend off;
end

% Fit imaginary part.
[xData, yData] = prepareCurveData(x, psi_i);
[fitP{2}, gof(2)] = fit(xData, yData, ft, opts);

coeff = coeffvalues(fitP{1}) + 1i*coeffvalues(fitP{2});
A = abs(coeff(2:end));          

phi_0 = angle(coeff(1));
phi = zeros(1,ord);
for i=1:ord
    phi(i) = angle(coeff(i+1));
    A(i) = A(i)*exp(1i*(phi(i) - phi_0));
end

if approx == 1
    d = str2double(handles.decimalsEdit.String);
    Dr = d - ceil(log10(abs(real(A))));
    Di = d - ceil(log10(abs(imag(A))));
    D = zeros(1, ord);
    for i=1:ord
        if Dr > Di
            D(i) = Di(i);
        else
            D(i) = Dr(i);
        end
    end
    Ar = zeros(1, ord); Ai = zeros(1, ord);
    for i=1:ord
        Ar(i) = round(real(A(i)), D(i));
        Ai(i) = round(imag(A(i)), D(i));
    end
    A = Ar + 1i*Ai;
end   

A02 = 1;
for i = 1:ord
    A02 = A02 - 2*abs(A(i))^2;
end
A0 = sqrt(A02);

% Plot fit with data.
if plot_res
    % plot fit
    subplot(2, 2, 3, 'Parent', handle);
    h = plot(fitP{2}, xData, yData );
    %legend(h, 'Data', 'Fit', 'Orientation', 'Horizontal', 'Location', 'Best');
    legend off;
    xlabel('x'); ylabel('\psi_r'); grid on;
    % plot residuals
    subplot(2, 2, 4, 'Parent', handle);
    h = plot(fitP{2}, xData, yData, 'Residuals');
    xlabel('x'); ylabel('\psi_r'); grid on;
    %legend(h, 'Residuals', 'Zero-line', 'Orientation', 'Horizontal', 'Location', 'Best');
    legend off;
end