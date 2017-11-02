function res=generatePdfLogNormal(S0,m,s,T,x)
m=(m-0.5*s^2)*T+log(S0)
s=s*sqrt(T)
res=exp(-((-m + log(x)).^2/(2*s^2)))./(sqrt(2*pi)*s*x);
end