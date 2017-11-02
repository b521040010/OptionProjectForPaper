function testInterestRatesOptimization()
    delete(findall(0,'Type','figure'));
    initutilityoptimization();
    riskAversion = 0.001;
    utilityFunction = ExponentialUtilityFunction( riskAversion );
    date = '20160408T150000';
    dayData = DayData( date );
    
    %same interest rates for lending and borrowing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %lending
    dayData.instruments{1}.r=0;
    %borrowing
    dayData.instruments{2}.r=0.001;
    %arbitrage = ArbitrageFinder.findArbitrageForDate( date, false, true );
    %assert(~arbitrage);
    model = dayData.blackScholesModel();
    model.mu = 0.08;
    
    % Create a utility maximization problem corresponding
    % to this problem
    ump = UtilityMaximizationProblem1D();
    ump.setModel( model );
    ump.setUtilityFunction(utilityFunction);
    zcb = dayData.findInstrument(0, DayData.cashType );
    %my own adjustment from %ump.setCurrentPosition(zcb);
    currentPort=Portfolio();
    currentPort.add([100000],{zcb})
    ump.setCurrentPosition(currentPort);
    %------------------------
    
    for i=1:length(dayData.instruments)
        ump.addInstrument( dayData.instruments{i} );
    end

    for idx = 1:length(ump.instruments)
        instrument=ump.instruments{idx};
        ump.addConstraint(QuantityConstraint(idx,-instrument.bidSize,instrument.askSize));
    end
    
    ump.addConstraint( BoundedLiabilityConstraint());
    [utility, quantities] = ump.optimize();  
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %lending
    dayData.instruments{1}.r=0;
    %borrowing
    dayData.instruments{2}.r=0.1;
    %arbitrage = ArbitrageFinder.findArbitrageForDate( date, false, true );
    %assert(~arbitrage);
    model = dayData.blackScholesModel();
    model.mu = 0.08;
    
    % Create a utility maximization problem corresponding
    % to this problem
    ump = UtilityMaximizationProblem1D();
    ump.setModel( model );
    ump.setUtilityFunction(utilityFunction);
    zcb = dayData.findInstrument(0, DayData.cashType );
    %my own adjustment from %ump.setCurrentPosition(zcb);
    currentPort=Portfolio();
    currentPort.add([100000],{zcb})
    ump.setCurrentPosition(currentPort);
    %------------------------
    
    for i=1:length(dayData.instruments)
        ump.addInstrument( dayData.instruments{i} );
    end

    for idx = 1:length(ump.instruments)
        instrument=ump.instruments{idx};
        ump.addConstraint(QuantityConstraint(idx,-instrument.bidSize,instrument.askSize));
    end
    
    ump.addConstraint( BoundedLiabilityConstraint());
    [utility2, quantities2] = ump.optimize();    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %lower rate for borrowing must be better for investors
    assert(utility>utility2)
    %higher interest rate for borrowing will discurage the borrowing
    assert(quantities(2)<quantities2(2))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %If the lending rate is higher than the borrowing rate, an arbitrage
    %exists
    %lending
    dayData.instruments{1}.r=0.1;
    %borrowing
    dayData.instruments{2}.r=0.01;
    %arbitrage = ArbitrageFinder.findArbitrageForDate( date, false, true );
    %assert(~arbitrage);
    model = dayData.blackScholesModel();
    model.mu = 0.08;
    
    % Create a utility maximization problem corresponding
    % to this problem
    ump = UtilityMaximizationProblem1D();
    ump.setModel( model );
    ump.setUtilityFunction(utilityFunction);
    zcb = dayData.findInstrument(0, DayData.cashType );
    %my own adjustment from %ump.setCurrentPosition(zcb);
    currentPort=Portfolio();
    currentPort.add([100000],{zcb})
    ump.setCurrentPosition(currentPort);
    %------------------------
    
    for i=1:length(dayData.instruments)
        ump.addInstrument( dayData.instruments{i} );
    end

    for idx = 1:length(ump.instruments)
        instrument=ump.instruments{idx};
        ump.addConstraint(QuantityConstraint(idx,-instrument.bidSize,instrument.askSize));
    end
    
    ump.addConstraint( BoundedLiabilityConstraint());
    [utility3, quantities3] = ump.optimize();  
    ump.plotPortfolio( sprintf('Net profit of the portfolio invested on %s. Utility=%d', date,utility),quantities3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

end