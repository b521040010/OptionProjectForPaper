function testInstruments()
%TESTINSTRUMENTS Summary of this function goes here
%   Detailed explanation goes here

bond = Bond(0.5,0.1,1,1,Inf,5000);
assert(bond.getBid()==1);
assert(bond.getAsk()==1);
assert(bond.getBidSize()==Inf);
assert(bond.getAskSize()==5000);
assertApproxEqual(bond.payoff(1),exp(0.5*0.1)*1, 1e-6);
newBond = bond.updateMaturity(1);
assertApproxEqual( newBond.T, 0.5 - 1/365, 0.0000001);
assertApproxEqual( newBond.getBid(), exp(0.1*1/365)*bond.getBid(), 0.0000001);
assertApproxEqual( newBond.getAsk(), exp(0.1*1/365)*bond.getAsk(), 0.0000001);

future = Future2(99,0,0,100,50);
assertApproxEqual(future.getBid(),0,1e-6);
assertApproxEqual(future.getAsk(),0,1e-6);
assertApproxEqual(future.payoff(105),6*50, 1e-6);
assertApproxEqual(future.payoff(50),-49*50, 1e-6);
assertApproxEqual(future.getBidSize(),100,1e-6);
assertApproxEqual(future.getAskSize(),50,1e-6);

% future = Future(99,101);
% assertApproxEqual(future.getBid(),99,1e-6);
% assertApproxEqual(future.getAsk(),101,1e-6);
% assertApproxEqual(future.payoff(105),105, 1e-6);
% future.updateMaturity(1);

call = CallOption(110, 1,2,20,30);
assertApproxEqual(call.payoff(105),0, 1e-6);
assertApproxEqual(call.payoff(115),5*100, 1e-6);
assertApproxEqual(call.getBid(),1*100-call.commission, 1e-6);
assertApproxEqual(call.getAsk(),2*100+call.commission, 1e-6);
assertApproxEqual(call.getBidSize(),20, 1e-6);
assertApproxEqual(call.getAskSize(),30, 1e-6);
newCall = call.updateMaturity(1);
assertApproxEqual( newCall.K, call.K, 0.0001);

put = PutOption(110, 1,2,Inf,60);
assertApproxEqual(put.payoff(105),5*100, 1e-6);
assertApproxEqual(put.payoff(115),0, 1e-6);
assertApproxEqual(put.getBid(),1*100-put.commission, 1e-6);
assertApproxEqual(put.getAsk(),2*100+put.commission, 1e-6);
assert(put.getBidSize()==Inf);
assertApproxEqual(put.getAskSize(),60, 1e-6);
newPut = put.updateMaturity(1);
assertApproxEqual( newPut.K, put.K, 0.0001 );


end

