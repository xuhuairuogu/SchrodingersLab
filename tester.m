t = linspace(0, 43.35, 100000);
r = sqrt(3)/2;
c = -38;
k = [0:255];
[tg, kg] = meshgrid(t, k);
A_k = exp(abs(tg).*kg.*r+kg*c);
psi = sum(A_k);
A_0 = sqrt(1 - (sum(A_k.^2)));
%figure
%plot(t, log(A_k), '*');

a = 3/8;
b = sqrt(8*a*(1-2*a));
omega = 2*sqrt(1-2*a);
%t_2 = t+43.7120;
t_2 = t-43.35;
x_2 = 0;
psi_2 = (1 + (2*(1-2*a)*cosh(b*t_2) + 1i*b*sinh(b*t_2))./(sqrt(2*a)*cos(omega*x_2)-cosh(b*t_2))).*exp(1i*t_2);

plot(t, abs(psi).^2, '-b', 'LineWidth', 2);
hold on;
plot(t, abs(psi_2).^2, '-r', 'LineWidth', 2);
legend('Reproduced breather', 'Analytical breather', 0);
xlim([35, max(t)]);
ylabel('|\psi|^2'); xlabel('t');
