function testInterestRatesOptimization()
    % This test is designed to test if the different rates for borrowing 
    % and lending is handled properly
    riskAversion = 0.001;
    utilityFunction = ExponentialUtilityFunction( riskAversion );
    date = '20160408T150000';
    dayData = DayData( date );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First, we will use the idea that lower rate for borrowing must be 
    %better for investors and higher interest rate for borrowing will 
    %discurage the borrowing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %lending rate
    dayData.instruments{1}.r=0;
    %borrowing rate
    dayData.instruments{2}.r=0.001;
    model = dayData.blackScholesModel();
    model.mu = 0.08;
    
    % Create a utility maximization problem corresponding
    % to this problem
    ump = UtilityMaximizationProblem1D();
    ump.setModel( model );
    ump.setUtilityFunction(utilityFunction);
    %set initial wealth
    zcb = dayData.findInstrument(0, DayData.cashType );
    currentPort=Portfolio();
    currentPort.add([100000],{zcb})
    ump.setCurrentPosition(currentPort);  
    for i=1:length(dayData.instruments)
        ump.addInstrument( dayData.instruments{i} );
    end
    for idx = 1:length(ump.instruments)
        instrument=ump.instruments{idx};
        ump.addConstraint(QuantityConstraint(idx,-instrument.bidSize,instrument.askSize));
    end    
    ump.addConstraint( BoundedLiabilityConstraint());
    [utility, quantities] = ump.optimize();  
    
    %lending rate
    dayData.instruments{1}.r=0;
    %borrowing rate
    dayData.instruments{2}.r=0.1; % This is higher than the previous one
    model = dayData.blackScholesModel();
    model.mu = 0.08;
    
    % Create a utility maximization problem corresponding
    % to this problem
    ump = UtilityMaximizationProblem1D();
    ump.setModel( model );
    ump.setUtilityFunction(utilityFunction);
    zcb = dayData.findInstrument(0, DayData.cashType );
    currentPort=Portfolio();
    currentPort.add([100000],{zcb})
    ump.setCurrentPosition(currentPort);
  
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
    %lending rate
    dayData.instruments{1}.r=0.1; %lending rate is higher
    %borrowing rate
    dayData.instruments{2}.r=0.01; 
    model = dayData.blackScholesModel();
    model.mu = 0.08;
    
    % Create a utility maximization problem corresponding
    % to this problem
    ump = UtilityMaximizationProblem1D();
    ump.setModel( model );
    ump.setUtilityFunction(utilityFunction);
    zcb = dayData.findInstrument(0, DayData.cashType );
    currentPort=Portfolio();
    currentPort.add([100000],{zcb})
    ump.setCurrentPosition(currentPort);
    
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
    %we can observe the arbitrage from the graph
    assert(utility3>utility2)
    assert(utility3>utility)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % If we only have a bank account with interest rate for lending=0.01 
    %and interest rate for borrowing=0.05,
    % We must to trade and the utilityty must be
    % (1-exp(-0.001*10^5*exp(0.01)))/0.001 (=E[v(100,000)])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dayData.clearInstruments;
    %lending
    dayData.addInstrument( Bond( 1, 0.01, 1,1,0,Inf));
    %borrowing
    dayData.addInstrument( Bond( 1, 0.05, 1, 1,Inf,0 ));

    model = BlackScholesModel();
    model.T=1;
    model.sigma=0.35;
    model.mu = 0.08;

    ump = UtilityMaximizationProblem1D();
    ump.setModel( model );
    ump.setUtilityFunction(utilityFunction);
    zcb = dayData.findInstrument(0, DayData.cashType );
    currentPort=Portfolio();
    currentPort.add([100000],{zcb})
    ump.setCurrentPosition(currentPort);
    
    for i=1:length(dayData.instruments)
        ump.addInstrument( dayData.instruments{i} );
    end
    
    for idx = 1:length(ump.instruments)
        instrument=ump.instruments{idx};
        ump.addConstraint(QuantityConstraint(idx,-instrument.bidSize,instrument.askSize));
    end
    
    ump.addConstraint( BoundedLiabilityConstraint());
    [utility4, quantities4] = ump.optimize();
    assert(utility4==(1-exp(-0.001*10^5*exp(0.01)))/0.001);
    assert(quantities4(1)==0);
    assert(quantities4(2)==0);
   
end