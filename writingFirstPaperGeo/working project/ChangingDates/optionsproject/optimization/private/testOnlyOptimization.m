function testOnlyOptimization()
    delete(findall(0,'Type','figure'));
    initutilityoptimization();
    riskAversion = 0.001;
    utilityFunction = ExponentialUtilityFunction( riskAversion );
    date = '20160408T145500';
    dayData = DayData( date );
    model = dayData.blackScholesModel();
    model.mu = 0.08;
    arbitrage = ArbitrageFinder.findArbitrageForDate( date, false, true );
    assert(~arbitrage);    
    % Create a utility maximization problem corresponding
    % to this problem
    ump = UtilityMaximizationProblem1D();
    ump.setModel( model );
    ump.setUtilityFunction(utilityFunction);
    % Set initial wealth
    zcb = dayData.findInstrument(0, DayData.cashType );
    currentPort=Portfolio();
    currentPort.add([100000],{zcb})
    ump.setCurrentPosition(currentPort);
    
    for i=1:length(dayData.instruments)
        ump.addInstrument( dayData.instruments{i} );
    end
    
    %add constraints
    for idx = 1:length(ump.instruments)
        instrument=ump.instruments{idx};
        ump.addConstraint(QuantityConstraint(idx,-instrument.bidSize,instrument.askSize));
    end    
    ump.addConstraint( BoundedLiabilityConstraint());
    
    [utility, quantities] = ump.optimize();
    ump.plotPortfolio( sprintf('Net profit of the portfolio invested on %s. Utility=%d', date,utility),quantities);
    plotPortfolio( ump.getInstruments(), quantities);
    
    %Test if buying and selling quantities are in [-bidSizes,askSizes]
    for idx=1:length(ump.instruments)
        instrument=ump.instruments{idx};
        assert(quantities(idx)<=instrument.askSize);
        assert(quantities(idx)>=-instrument.bidSize);
    end
    totalInvestment=0;
    for i = 1:length(ump.instruments)
        if quantities(i)>=0
            totalInvestment = totalInvestment+quantities(i)*ump.instruments{i}.getAsk();
        else
            totalInvestment = totalInvestment+quantities(i)*ump.instruments{i}.getBid();
        end
    end
    temp=values(ump.currentPosition.map);
    initialInvestment=temp{1}.quantity;
    % Test that the initial invested money is less than our initial wealth
    assert(totalInvestment<=initialInvestment);
end