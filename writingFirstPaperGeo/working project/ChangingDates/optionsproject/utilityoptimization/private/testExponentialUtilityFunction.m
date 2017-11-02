function testExponentialUtilityFunction()

eufRisky = ExponentialUtilityFunction(0.7);
eufRiskAverse = ExponentialUtilityFunction(1.5);
eufRisky.testEvaluation( 1.8 );
eufRiskAverse.testEvaluation( 1.8 );
assert( eufRisky.weightedUtility( 1.8, 0 ) > eufRiskAverse.weightedUtility( 1.8, 0 ) );

end

