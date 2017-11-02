function testDeltaHedgeProblem()
%TESTDELTAHEDGEPROBLEM Summary of this function goes here
%   Detailed explanation goes here

prob = DeltaHedgeProblem();
prob.proportionTransactionCosts = 0;
prob.bsm = BlackScholesModel();
prob.bsm.S0 = 1;
prob.K = 1;
prob.nPaths = 100;
price = prob.bsm.price( prob.r, prob.isPut, prob.K );
rng('default');
finalBalance = prob.simulateDeltaHedge( 1, price, prob.simulatePricePaths(1000) );
assert( std(finalBalance)>0.0 );
assertApproxEqual(mean(finalBalance),0, 0.1);

quantity = 1/prob.bsm.S0;
eu10 = prob.computeExpectedUtility( quantity, price*quantity, prob.simulatePricePaths(10) );
eu100 = prob.computeExpectedUtility( quantity, price*quantity, prob.simulatePricePaths(100) );
eu1000 = prob.computeExpectedUtility( quantity,price*quantity, prob.simulatePricePaths(1000) );
assert( eu1000<0.0 );
assert( eu100<eu1000 );
assert( eu10<eu100 );

nSteps = 10;
S = prob.simulatePricePaths(nSteps);
sellersPrice =prob.sellersIndifferencePrice( quantity, price*quantity, S );
assert( sellersPrice > price*quantity );
euWithSellersPrice = prob.computeExpectedUtility( quantity, sellersPrice, S );
assertApproxEqual( euWithSellersPrice, 0, 0.00001);

end

