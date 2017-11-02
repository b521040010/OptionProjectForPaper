function res=generatePdfLogST(S0,m,sigma,nu,x)
m=m+log(S0)
s = log(x);
            res = gamma((nu+1)/2)/(gamma(nu/2)*sqrt(pi*nu)*sigma) ...
                *( 1 + 1/nu*((s-m)/sigma).^2).^(-(nu+1)/2) ...
                ./ x;
end