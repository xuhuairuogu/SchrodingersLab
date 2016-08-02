function [fitresult, gof] = genCoeff(x, psi, a, ord, plot_res, handle)
%CREATEFITS(X_DT,PSI_R,PSI_I)
%  Create fits for function supplied by DT to be used as initial condition
%  for RW/AB of order "ord"
%
%  Input:
%       x: x-coordinated used with psi below.
%       psi: array containing complex psi from DT
%       a: a parameter from DT.
%       ord: order of RW, and thus, order of Fourier Series fit
%       plot_res: produce result of fit + residuals
%  Output:
%      fitresult : a cell-array of fit objects representing the fits.
%      gof : structure array with goodness-of fit info.

% Initialization.
Omega = 2*sqrt(1-2*a);
psi_i = imag(psi.'); psi_r = real(psi.');

% Initialize arrays to store fits and goodness-of-fit.
fitresult = cell(2, 1);
gof = struct( 'sse', cell(2,1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

% Set up fit type and options.
fitString = sprintf('A0 + 2*A1*cos(%.16f*x)', Omega); % First order fit
for i=2:ord % Higher order fits
    fitString = strcat(fitString, sprintf(' + 2*A%d*cos(%d*%.16f*x)', i, i, Omega));
end
ft = fittype(fitString, 'independent', 'x', 'dependent', 'y' );
opts = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', 0.5*ones(1,ord+1));
opts.Display = 'Off';

% Fit real part
[xData, yData] = prepareCurveData(x, psi_r);
[fitresult{1}, gof(1)] = fit(xData, yData, ft, opts);

% Plot fit with data.
if plot_res
    %figure('Name', 'Fits');
    % plot fit
    subplot(2, 2, 1, 'Parent', handle);
    h = plot(fitresult{1}, xData, yData );
    %legend(h, 'Data', 'Fit', 'Orientation', 'Horizontal', 'Location', 'Best');
    legend off;
    xlabel('x'); ylabel('\psi_r'); grid on;
    % plot residuals
    subplot(2, 2, 2, 'Parent', handle);
    h = plot(fitresult{1}, xData, yData, 'Residuals');
    xlabel('x'); ylabel('\psi_r'); grid on;
    %legend(h, 'Residuals', 'Zero-line', 'Orientation', 'Horizontal', 'Location', 'Best');
    legend off;
end

% Fit imaginary part.
[xData, yData] = prepareCurveData(x, psi_i);
[fitresult{2}, gof(2)] = fit(xData, yData, ft, opts);

% Plot fit with data.
if plot_res
    % plot fit
    subplot(2, 2, 3, 'Parent', handle);
    h = plot(fitresult{2}, xData, yData );
    %legend(h, 'Data', 'Fit', 'Orientation', 'Horizontal', 'Location', 'Best');
    legend off;
    xlabel('x'); ylabel('\psi_r'); grid on;
    % plot residuals
    subplot(2, 2, 4, 'Parent', handle);
    h = plot(fitresult{2}, xData, yData, 'Residuals');
    xlabel('x'); ylabel('\psi_r'); grid on;
    %legend(h, 'Residuals', 'Zero-line', 'Orientation', 'Horizontal', 'Location', 'Best');
    legend off;
end