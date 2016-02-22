function [psi] = AB_je(x, z, k)
    if k >= 1;
    	[SN, CN, DN] = ellipj(k*z, 1/k);
    	[~, CNA, ~] = ellipj(sqrt(2*k)*x, sqrt((k-1)/(2*k)));
    
    	A = sqrt(1/(1+k))*CNA;
    
    	psi = k*(A.*DN+1i/k*SN)./(1-A.*CN).*exp(1i*z);
    elseif k < 1
        [SN, CN, DN] = ellipj(z, k);
        [~, CNA, DNA] = ellipj(x*sqrt(1+k), sqrt((1-k)/(1+k)));
        
        A = CNA./DNA;
        
        psi = k*(A.*CN+1i*sqrt(1+k)*SN)./(sqrt(1+k)-A.*DN).*exp(1i*z);
   end
end
