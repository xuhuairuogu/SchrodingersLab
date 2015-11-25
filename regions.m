function [suprema] = regions(PSI, x, t)
[row,col] = find(abs(PSI).^2' >= 7);
results = zeros(length(row)+1, 3);
for i = 1:length(row)
    results(i, :) = [t(col(i)) x(row(i)) abs(PSI(col(i), row(i))).^2];
end
results(length(row)+1, :) = [9999999999999 0 0];

pivot = 1;
suprema = [];
for i = 1:length(row)+1;
    if results(i, 1) - results(pivot, 1) < 2.5
        continue
    else
        region = results(pivot:i-1, :);
        [~, index] = max(region(:, 3));
        suprema = [suprema; region(index, :)];
        pivot = i;
    end
end
end

