function [psi] = AB_je(x, z, k)

	[SN, CN, DN] = ellipj(k*z, 1/k);
	[~, CNA, ~] = ellipj(sqrt(2*k)*x, sqrt((k-1)/(2*k)));

	A = sqrt(1/(1+k))*CNA;

	psi = (k*A.*DN+1i/k*SN)./(1-A.*CN)*exp(1i*z);

end
