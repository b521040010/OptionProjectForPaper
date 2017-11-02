function testDayData()

dd = DayData( '20160408T150000' );

r = dd.getInterestRate();
T = dd.getT();
assertApproxEqual(T,70/365,1e-6);
assertApproxEqual( dd.bidPrice(2051.25,DayData.futureType),0,1e-6);
assertApproxEqual( dd.askPrice(2051.25,DayData.futureType),0,1e-6);
assertApproxEqual( dd.bidPrice(2051.5,DayData.futureType),0,1e-6);
assertApproxEqual( dd.askPrice(2051.5,DayData.futureType),0,1e-6);
callOptionStrike1200=dd.findInstrument( 1200, DayData.callType);
putOptionStrike1200=dd.findInstrument( 1200, DayData.putType);
assert( dd.bidPrice(1200,DayData.callType)==848.9*100-callOptionStrike1200.commission);
assert( dd.askPrice(1200,DayData.callType)==852.4*100+callOptionStrike1200.commission);
assert( dd.bidPrice(1200,DayData.putType)==0.05*100-putOptionStrike1200.commission);
assert( dd.askPrice(1200,DayData.putType)==0.6*100+putOptionStrike1200.commission);
assert(dd.bidSize(800,DayData.callType)==9);
assert(dd.askSize(800,DayData.callType)==11);
assert(dd.bidSize(1125,DayData.callType)==1);
assert(dd.askSize(1125,DayData.callType)==11);
assert(dd.bidSize(1825,DayData.callType)==10);
assert(dd.askSize(1825,DayData.callType)==9);
assert(dd.bidSize(2135,DayData.callType)==568);
assert(dd.askSize(2135,DayData.callType)==653);
assert(dd.bidSize(1825,DayData.putType)==607);
assert(dd.askSize(1825,DayData.putType)==610);
assert(dd.bidSize(2135,DayData.putType)==9);
assert(dd.askSize(2135,DayData.putType)==11);
assert(dd.bidSize(2051.25,DayData.futureType)==195);
assert(dd.askSize(2051.25,DayData.futureType)==0);
assert(dd.bidSize(2051.5,DayData.futureType)==0);
assert(dd.askSize(2051.5,DayData.futureType)==136);



% dd = DayData( '2014-01-16' );
dd = DayData('20160408T150000');
m = dd.blackScholesModel();
expectedT = datenum([2016 6 17 0 0 0])-datenum([2016 4 8 0 0 0]);
assertApproxEqual( m.T, expectedT/365.0, 1e-6);
assertApproxEqual( m.S0, 2058.65,0.0001);
assertApproxEqual(m.sigma,0.1329838,0.0001);
m = dd.calibrateHistoric( BlackScholesModel(), 0.90 );
assertApproxEqual(m.sigma,0.1174,0.001);
end

