function testIndifferencePricing()

    riskAversion = 0.02;
    utilityFunction = ExponentialUtilityFunction( riskAversion );
    date = '20160408T145920';
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
    %my own adjustment from %ump.setCurrentPosition(zcb);
    currentPort=Portfolio();
    currentPort.add([1000000],{zcb})
    ump.setCurrentPosition(currentPort);
    %------------------------
    

    for i=1:length(dayData.instruments)
        ump.addInstrument( dayData.instruments{i} );
    end
    
    ump.addConstraint(QuantityConstraint(2,-Inf,0));
    ump.addConstraint(QuantityConstraint(3,0,Inf));
    for idx = 1:length(ump.instruments)
        instrument=ump.instruments{idx};
        ump.addConstraint(QuantityConstraint(idx,-instrument.bidSize,instrument.askSize));
    end
    


    ump.addConstraint( BoundedLiabilityConstraint());
    [utility, quantities] = ump.optimize()
    ump.plotPortfolio( sprintf('Portfolio with only bounded liability constraint. Utility=%d', utility),quantities);
    plotPortfolio( ump.getInstruments(), quantities);
    n=length(quantities);
    
    %Test if buying and selling quantities are in [-bidSizes,askSizes]
    for idx=1:length(ump.instruments)
        instrument=ump.instruments{idx};
        assert(quantities(idx)<=instrument.askSize);
        assert(quantities(idx)>=-instrument.bidSize);
    end
    
    % This means that we can at least solve the utility maximization
    % problem to begin with
% 
%     K = 1940;
%     call = dayData.findInstrument( K, DayData.callType );
%     put = dayData.findInstrument( K, DayData.putType );
%     ump.removeInstrument( call );
%     ump.removeInstrument( put );
%     
%     % Get rid of instruments which are too far out of the money
%     
% %      function accept = filter( instrument )
% %          if (isa(instrument,'CallOption') || isa(instrument,'PutOption'))
% %              accept = abs( instrument.K - K) < 50;
% %          else 
% %              accept = 1;
% %          end
% %      end
% %      ump.filterInstruments( @filter );
% %     
%     proportions = 0.5:-0.1:-0.2;
%     indifferencePrices = zeros( 1, length( proportions ));
%     zcb = dayData.findInstrument(0, DayData.cashType );
%     for i=1:length( proportions)
%         price = call.getAsk();
%         proportion = proportions(i);
%         try 
%             indifferencePrices(i) = ump.indifferencePrice(zcb, proportion/price, call, proportion);        
%         catch ex
%             indifferencePrices(i) = NaN;
%         end
%     end
% 
%     figure();
%     plot( proportions, indifferencePrices );
%     xlabel('Proportion of $1 ask quantity');
%     ylabel('Indifference price');
%     title( sprintf('Indifference price of call option with strike %d', K));
end

