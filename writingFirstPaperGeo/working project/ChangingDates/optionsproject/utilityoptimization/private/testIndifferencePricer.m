function testIndifferencePricer()

S0 = 100;
meanT = 108;
sdT = 20;
T = 1;
r = 0;
gamma = 3;

bm = BachelierModel();
bm.S0 = S0;
bm.meanT = meanT;
bm.sdT = sdT;

prob = UtilityMaximizationProblem1D();
prob.setModel(bm);
zcb = Bond(T,r,1,1,Inf,Inf);

%%%%%
currentPort=Portfolio();
currentPort.add([0],{zcb});
prob.setCurrentPosition(currentPort);
%%%%%%%%

%%%%%%%%%%%%
% currentPort=Portfolio();
% currentPort.add([0],{zcb})
% prob.addInstrument( currentPort );
%%%%%%%%%%%%%%%%%%
prob.addInstrument( zcb );
prob.setUtilityFunction( ExponentialUtilityFunction( gamma ));

instrument = Future(exp(r*T)*S0,exp(r*T)*S0,Inf,Inf);
%%%%%%%%%%%%%%
bond=Bond(T, r, 1,1,0,0);
%%%%%%%%%%%%%
% If the instrument is in the market then the only possible price
% is the instrument price
copyProb = prob.copy();
%%%%%%%%%%
% copyProb.addInstrument(bond);
%%%%%%%%%%%
copyProb.addInstrument(instrument);
price = copyProb.indifferencePrice( zcb, 1/S0, instrument, 1 );
assertApproxEqual(price,1, 1e-4);
copyProb.addConstraint( QuantityConstraint(2,-2,2));
price = copyProb.indifferencePrice( zcb, -1/S0, instrument, 1 );
assertApproxEqual(price,-1, 1e-4);

p1 = prob.indifferencePrice( zcb, 0.5*1/S0, instrument, 1 );
p2 = prob.indifferencePriceExponentialUtility( zcb, 0.5*1/S0, instrument, 1 );
assertApproxEqual(p1,p2,1e-4);

proportions = -5:0.5:5;
prices = zeros(1,length( proportions));
for i=1:length(proportions)
    prices(i) = prob.indifferencePrice( zcb, proportions(i)*1/S0, instrument, 1 );
end
figure();
plot( proportions, prices );
title('Indifference price of stock - Bachelier Model');
xlabel('Market price');
ylabel('Indifference price');

% Test the no arbitrage bounds
[sup,sub] = prob.superHedgePrice( zcb );
assertApproxEqual( sup, 1.0, 0.0001 ); 
assertApproxEqual( sub, 1.0, 0.0001 ); 

end

