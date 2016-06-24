%close all
clear variables

a0 = 0.43;
da = 0.01;
af = 0.43;
ak = a0:da:af;
T = 15;                                 
xj = [0, 0, 0, 0, 0];
k0 = 1;
kf = 45;
height = zeros(1, kf);
base = zeros(1, kf);
area = zeros(1, kf);
for i = 1:length(ak);
    a = ak(i);
    for j = k0:kf
        tj = [-j/5-1, 0, 0, 0, 0];
        [PSI, x, t] = calcDarboux(2, a, 1024, 1000, T, xj, tj);
        maxima = regions(PSI, x, t);
        height(j) = maxima(2,1) - maxima(1,1);
        base(j) = abs(maxima(2,2))*2;
        area(j) = 0.5*base(j)*height(j);
    end
    figure;
    Omega = 2*sqrt(1-2*ak);
    L = pi/sqrt(1-2*a);
    hold all;
    plot([1, 10], [L/2, L/2], '-k', 'LineWidth', 2);
    l = 6:0.2:10;
    h = height(25:end);
    [fitresult, ~] = fit(l.', h.', 'poly1');
    coeff = coeffvalues(fitresult); c = coeff(2);
    l2 = 1.1:0.1:10;
    y = l2 + c;
    plot(l2, y, '-b', 'LineWidth', 2);
    string = sprintf('t_{s1} + %.3f', c);
    plot((k0:kf)/5+1, height, 'rx', 'MarkerSize', 8); grid on;
    plot((k0:kf)/5+1, base, '+', 'MarkerEdgeColor', [0, 0.4980, 0]);
    hold off;
    legend('L/2', string, 'Base', 'Height', 'Location', 'Northwest'); 
    title(sprintf('a = %.3f', a)); xlabel('t_{s1}'); ylabel('Spacing');
    ylim([2, 14]);
end