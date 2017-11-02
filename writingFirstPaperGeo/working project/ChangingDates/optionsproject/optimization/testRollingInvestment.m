function testRollingInvestment()
%TESTROLLINGINVESTMENT Summary of this function goes here
%   Detailed explanation goes here
initutilityoptimization();

startDate = '2014-01-02';
endDate = '2014-10-21';

ri = RollingInvestment( startDate, endDate, 0.0);
ri.run();



end

