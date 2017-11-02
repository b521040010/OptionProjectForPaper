try
    %riskAversion=[0.000001:0.000001:0.00003];
    %riskAversion=[0.00003:-0.000001:0.000001];
    riskAversion=0.00002;
    %K=[2000:50:3000];
    quantity=1;
    sigma=[0.1:0.0025:0.2];
    K=2000;
    n=length(sigma);
    indifferencePricesForBuying=zeros(1,n);
    indifferencePricesForSelling=zeros(1,n);
    sup=zeros(1,n);
    sub=zeros(1,n);
    for i=1:n
        sigma(i)
        [indifferencePricesForBuying(i),~,~,~]=testComputeIndifferencePrice(riskAversion,K,quantity,sigma(i));
        [indifferencePricesForSelling(i),~,~,~]=testComputeIndifferencePrice(riskAversion,K,-quantity,sigma(i));
        indifferencePricesForBuying
        indifferencePricesForSelling
        %sup
        %sub
    end
    indifferencePricesForBuying=indifferencePricesForBuying./(quantity*100);
    indifferencePricesForSelling=indifferencePricesForSelling./(-quantity*100);

         
catch
  x=1;

end
    