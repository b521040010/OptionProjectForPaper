function testMarketData()
%TESTMARKETDATA Summary of this function goes here
%   Detailed explanation goes here

marketData = MarketData.getInstance();
date = datenum( '08/04/2016', 'dd/mm/yyyy');
interestRate = marketData.getInterestRate(date);
assertApproxEqual(interestRate, log(1.004347), 0.0000001);
date = datenum( '12/03/2016', 'dd/mm/yyyy');
interestRate = marketData.getInterestRate(date);
assertApproxEqual(interestRate, log(1 + 2/3 * 0.004362 + 1/3 * 0.004413), 0.000000001);

set = marketData.getSPXSettlementValue( datenum('2016-02-19'));
assertApproxEqual(set, 1911.37, 0.000000001);


end

