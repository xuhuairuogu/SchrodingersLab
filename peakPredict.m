function peak = peakPredict(a, k)

%peak = -(1 - 4*a + sqrt(2*a))./(sqrt(2*a) - 1);
peak = -1;
for i=0:k
    peak = peak + 2*sqrt(2*i^2*a - i^2 + 1);
end

end
