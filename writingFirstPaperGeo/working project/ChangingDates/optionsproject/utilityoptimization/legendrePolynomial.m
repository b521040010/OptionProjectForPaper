function coeffs = legendrePolynomial( m )
%LEGENDREPOLYNOMIAL Generates the given legendre polynomial
if m==0
    coeffs = 1;
elseif m==1
    coeffs = [0 1];
else
    n = m-1;
    pn = legendrePolynomial(m-1);
    pnm1 = legendrePolynomial(m-2);
    coeffs = 1/(n+1)*((2*n+1)*[0 pn] - n*[pnm1 0 0]);
end        

end

