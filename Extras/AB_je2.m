function [psi] = AB_je2(x, z, k)

    k = 1/k;
    x = sqrt(2)*x;
	[SN, CN, DN] = ellipj(z/k, k);
	[~, CNA, ~] = ellipj(x/sqrt(k), sqrt((1-k)/2));

	A = sqrt(k/(1+k))*CNA;

	psi = (A.*DN+1i*k*SN)./(1-A.*CN).*exp(1i*z)/sqrt(2)/k;
    psi = sqrt(2)*psi;
end
