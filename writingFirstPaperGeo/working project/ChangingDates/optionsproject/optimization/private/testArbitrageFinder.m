function testArbitrageFinder()
%TESTARBITRAGEFINDER Summary of this function goes here
%   Detailed explanation goes here
%instruments=dayData.instruments; allowBonds=0; allowFutures=1;
global plotsInTests

model = BlackScholesModel();
model.sigma = 0.1;

dayData = DayData( '20160408T150000' );
dayData.clearInstruments();
dayData.spot=100;
dayData.addInstrument( Bond(1,0,1,1,Inf,Inf));
dayData.addInstrument( Future2(100,0,0,Inf,0));
dayData.addInstrument( Future2(100,0,0,0,Inf));
dayData.addInstrument( CallOption(100,20,20,Inf,Inf));
dayData.addInstrument( PutOption(100,20,20,Inf,Inf));

arbitrageFinder = ArbitrageFinder(dayData);
[arbitrage ] = arbitrageFinder.findArbitrage(false,true);
assert(~arbitrage);

% Bid ask arbitrage

dayData.clearInstruments();
dayData.addInstrument( Bond(1,0,1,1,Inf,Inf));
%dayData.addInstrument( Future(99,101,Inf,Inf));
 dayData.addInstrument( Future2(99,0,0,Inf,0));
 dayData.addInstrument( Future2(101,0,0,0,Inf));
dayData.addInstrument( CallOption(100,21,20,Inf,Inf));
dayData.addInstrument( PutOption(100,19,21,Inf,Inf));

arbitrageFinder = ArbitrageFinder(dayData);
[arbitrage ] = arbitrageFinder.findArbitrage(false,true);
assert(arbitrage);

% % Put call parity arbitrage 1
dayData.clearInstruments();
dayData.addInstrument( Bond(1,0,1,1,Inf,Inf));
dayData.addInstrument( Future2(105,0,0,Inf,0));
dayData.addInstrument( Future2(105,0,0,0,Inf));
dayData.addInstrument( CallOption(100,20,20,Inf,Inf));
dayData.addInstrument( PutOption(100,20,20,Inf,Inf));
arbitrageFinder = ArbitrageFinder(dayData);
[arbitrage ] = arbitrageFinder.findArbitrage(true,true);
assert(arbitrage);



end

