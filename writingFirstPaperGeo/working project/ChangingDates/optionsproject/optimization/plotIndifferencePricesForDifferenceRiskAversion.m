function [indifferencePricesForBuying]=plotIndifferencePricesForDifferenceRiskAversion()
riskAversion=[0.00006:0.000001:0.0002];
K=2000;
quantity=-1;
n=length(riskAversion);

indifferencePricesForBuying=zeros(1,n);
indifferencePricesForSelling=zeros(1,n);
%sup=zeros(1,n);
%sub=zeros(1,n);
for i=1:n
%[indifferencePricesForBuying(i),~,sup(i),sub(i)]=testComputeIndifferencePrice(riskAversion(i),K,quantity);
riskAversion(i)
[indifferencePricesForBuying(i),~]=testComputeIndifferencePrice(riskAversion(i),K,quantity);
end
indifferencePricesForBuying=indifferencePricesForBuying./(quantity*100);
%sup=sup./(quantity*100);
%sub=sub./(quantity*100);

% for i=1:n
% [indifferencePricesForSelling(i),~,~,~]=testComputeIndifferencePrice(riskAversion(i),K,-quantity);
% [indifferencePricesForSelling(i),~]=testComputeIndifferencePrice(riskAversion(i),K,-quantity);
% 
% end
% indifferencePricesForSelling=indifferencePricesForSelling./(-quantity*100);

%optional
%plot(K,-indifferencePrices/100)
%plot the price from Bloomberg of the asset we are pricing 
% bloombergPrice=
% arrayBloombergPrice=bloombergPrice*ones(1,n)
end
