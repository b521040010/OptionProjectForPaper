function testIndifferencePricingWithBidAskQuantities()
    delete(findall(0,'Type','figure'));
    initutilityoptimization();
    riskAversion = 2;
    utilityFunction = ExponentialUtilityFunction( riskAversion );
    date = '2014-10-19';

    dayData = DayData( date );
    arbitrage = ArbitrageFinder.findArbitrageForDate( date, false, true );
    %assert(~arbitrage);

    model = dayData.blackScholesModel();
    model.mu = 0.08;
    
    % Create a utility maximization problem corresponding
    % to this problem
    ump = UtilityMaximizationProblem1D();
    ump.setModel( model );
    ump.setUtilityFunction(utilityFunction);
    zcb = dayData.findInstrument(0, DayData.cashType );
    ump.setCurrentPosition(zcb);

    for i=1:length(dayData.instruments)
        ump.addInstrument( dayData.instruments{i} );
    end
    %{
    for idx=2:length(dayData.instruments)
        ump.addConstraint(QuantityConstraint(idx,-100,100));
    end
    %}
    ump.addConstraint( BoundedLiabilityConstraint());
    [utility, quantities] = ump.optimize();
    ump.plotPortfolio( sprintf('Portfolio with only bounded liability constraint. Utility=%d', utility),quantities);
    plotPortfolio( ump.getInstruments(), quantities);

    
    % This means that we can at least solve the utility maximization
    % problem to begin with

    K = 1940;
    K = 2060
     
    call = dayData.findInstrument( K, DayData.callType );
    put = dayData.findInstrument( K, DayData.putType );
    ump.removeInstrument( call );
    ump.removeInstrument( put );
    
    
    % Get rid of instruments which are too far out of the money
    function accept = filter( instrument )
        if (isa(instrument,'CallOption') || isa(instrument,'PutOption'))
            accept = abs( instrument.K - K) < 50;
        else 
            accept = 1;
        end
    end
    ump.filterInstruments( @filter );
    
    proportions = 0.5:-0.1:-0.2;
    
    indifferencePrices = zeros( 1, length( proportions ));
    zcb = dayData.findInstrument(0, DayData.cashType );
    for i=1:length( proportions)
        price = call.getAsk();
        proportion = proportions(i);

        try 
            indifferencePrices(i) = ump.indifferencePrice(zcb, proportion/price, call, proportion)      

        catch ex
            indifferencePrices(i) = NaN;
        end
    end

    figure();
    plot( proportions, indifferencePrices );
    xlabel('Proportion of $1 ask quantity');
    ylabel('Indifference price');
    title( sprintf('Indifference price of call option with strike %d', K));
    
    
end

