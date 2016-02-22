function [psi] = AB_je_new(x, z, k)

    x = sqrt(2)*x;
	[SN, CN, DN] = ellipj(z, k);
	[~, CNA, DNA] = ellipj(x*sqrt((1+k)/2), sqrt((1-k)/(1+k)));

	A = CNA./DNA;

	psi = k/sqrt(2)*(A.*CN+1i*sqrt(1+k)*SN)./(sqrt(1+k)-A.*DN).*exp(1i*z);
    psi = sqrt(2)*psi;
end
