function [suprema] = regions(PSI, x, t)
[peaks, imax] = extrema2(abs(PSI).^2);
[i, j] = ind2sub(size(PSI), imax);

t_pos = t(i);
x_pos = x(j);

suprema = [t_pos, x_pos, peaks];

suprema = sortrows(suprema);

end
