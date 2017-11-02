function testoptimization()

    delete(findall(0,'Type','figure'));
    initutilityoptimization();
     testArbitrageFinder();
     testDeltaHedgeProblem();
                % we dont use this class for now
                %testOptimizationProblem();
     testSPX();
     testDayData();
     testMarketData(); 
     testOnlyOptimization();
     testOnlyIndifferencePricing();
     testOldPortfolio(); 
     testExponentialUtility();
    %testExponentialOptimizationProblem();
     testInterestRatesOptimization();

end

