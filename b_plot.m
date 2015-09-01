function [A_0] = b_plot(PSI, Nx, t, num, peakTime)

J = 50;                                    % Scaling down of number of points
PSI_k = abs(fft(PSI'))/Nx;                  % Absolute normalized fft
data = log(PSI_k(1:num, 1:J:end));          % Our data

[tg, kg] = meshgrid(t, 1:num-1);

q = 3/8;
a = sqrt(8*q*(1-2*q));
alpha = sqrt(2*q)./cosh(a*(tg-peakTime));
alpha_2 = sqrt(2*q)./cosh(a*(t-peakTime));

A_0 = 1 - ((2*(1-2*q)+1i*a*tanh(a*(t-peakTime)))./sqrt(1-alpha_2.^2));
A_k = -(2*(1-2*q)+1i*a*tanh(a*(tg-peakTime)))./sqrt(1-alpha.^2).*((1-sqrt(1-alpha.^2))./alpha).^abs(kg);

plot(t, log(abs([A_0; A_k])), '-', 'LineWidth', 1.5);
hold on
plot(t(1:J:end), data, 'o', 'MarkerSize', 6); grid on;
ylim([-50, 5]);
xlim([0, max(t)]);
xlabel('t');
ylabel('log(|A_k|)');