function testPowerUtilityFunction()
%TESTPOWERUTILITYFUNCTION Summary of this function goes here
%   Detailed explanation goes here

pufRisky = PowerUtilityFunction(0.7);
pufLog = PowerUtilityFunction(1.0);
pufRiskAverse = PowerUtilityFunction(1.5);
pufRisky.testEvaluation( 1.8 );
pufLog.testEvaluation( 1.8 );
pufRiskAverse.testEvaluation( 1.8 );
assert( pufRisky.weightedUtility( 1.8, 0 ) > pufLog.weightedUtility( 1.8, 0 ) );
assert( pufLog.weightedUtility( 1.8, 0 ) > pufRiskAverse.weightedUtility( 1.8, 0 ) );


end

