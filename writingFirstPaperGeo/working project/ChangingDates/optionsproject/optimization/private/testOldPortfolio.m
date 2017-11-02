function testOldPortfolio()

p = OldPortfolio();
p.strikes = [100 110 120 100];
p.quantities = [10 20 30 -5];
put = DayData.putType;
call = DayData.callType;
p.instrumentTypes = [put call call put];

p = p.simplify();
assert( length( p.strikes)==3);
assert( sum( p.quantities )==55);
 
% Test payoff function
payoff = p.payoff(0.05, 115);
expectedPayoff = 20*(115-110);
assert(payoff==expectedPayoff);
% 
% Payoff of in the money put
p.strikes = 100;
p.quantities = -2;
p.instrumentTypes = put;
assert( p.payoff(0.05, 95) == -2*5);
% 
% Payoff of stock
p.strikes = 0;
p.quantities = 1;
p.instrumentTypes = DayData.futureType;
assert( p.payoff(0.05, 95) == 95);
p.quantities = -1;
assert( p.payoff(0.05, 95) == -95);
% 
% Payoff of cash
p.strikes = 0;
p.quantities = 1;
p.instrumentTypes = DayData.cashType;
assert( p.payoff(0.05, 95) == 1.05);
p.quantities = -1;
assert( p.payoff(0.05, 95) == -1.05);
% 
% 
% Test pricing functions
% 
dayData = DayData( '20160408T150000' );
% 
dayData.clearInstruments();
T = dayData.getT();
r = dayData.getInterestRate();
dayData.addInstrument( Bond( T,r,1,1,Inf,Inf));
dayData.addInstrument( Future2( 999.9,0,0,Inf,Inf));
dayData.addInstrument( CallOption( 1010, 20.2, 20.4, Inf, Inf ));
dayData.addInstrument( PutOption( 1010, 26.1, 26.2, Inf, Inf ));
% 
assertApproxEqual( dayData.askPrice(999.9,DayData.futureType),0,1e-6);
% 
csh = DayData.cashType;
future = DayData.futureType;
call = DayData.callType;
put = DayData.putType;
p.quantities = [500 1 2 -5];
p.strikes = [0 999.9 1010 1010];
p.instrumentTypes = [ csh future call put ];
callCommission=CallOption( 1010, 20.2, 20.4, Inf, Inf ).commission;
putCommission=PutOption( 1010, 26.1, 26.2, Inf, Inf ).commission;
assert( p.cost( dayData )==500 + 0 + 2*(20.4*100+callCommission) -5*26.1*100 +5*putCommission);
assert( p.value( dayData )==500 + 0 + 2*20.2*100+2*callCommission -5*26.2*100 +5*putCommission );

end

