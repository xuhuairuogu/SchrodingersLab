a1 = [0:0.001:0.5];
a2 = a1;
%[a1, a2] = meshgrid(a1, a2);
for i = 1:length(a1);
    for j = 1:length(a2);
        if a1(i) ~= a2(j)
            peak(i, j) = 1+2*sqrt(2*a1(i))+2*sqrt(2*a2(j));
        else
            peak(i,j) = 0;
        end
    end
end
surf(a1, a2, peak, 'EdgeColor', 'none');
figure; surf(a1, a2, peak, 'EdgeColor', 'none');

