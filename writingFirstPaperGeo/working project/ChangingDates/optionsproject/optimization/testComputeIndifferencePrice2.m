function [indifferencePrice,sup,sub,quantities]=testComputeIndifferencePrice2(riskAversion,K,quantity)
%[indifferencePrice,sub,sup,quantities]=testComputeIndifferencePrice(riskAversion,K,quantity)    
%delete(findall(0,'Type','figure'));

    initutilityoptimization();
    utilityFunction = ExponentialUtilityFunction( riskAversion );
    date = '20160408T145500';
    dayData = DayData( date );    
    %arbitrage = ArbitrageFinder.findArbitrageForDate( date, false, true );
    %assert(~arbitrage);

    %model = dayData.studentTModel();
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
    
    ump.addConstraint( BoundedLiabilityConstraint());

%      digitalCall=DigitalCallOption(K,1, 1.2 ,10,10)
    call = dayData.findInstrument( K, DayData.callType );
    put = dayData.findInstrument( K, DayData.putType );
    ump.removeInstrument( call );
    ump.removeInstrument( put );
     
%     Get rid of instruments which are too far out of the money
%     ump.instruments
%      function accept = filter( instrument )
%          if (isa(instrument,'CallOption') || isa(instrument,'PutOption'))
%              accept = abs( instrument.K - K) < 100;
%          else 
%              accept = 1;
%          end
%      end
%      ump.filterInstruments( @filter );
%     ump.instruments

     ump.addConstraint(QuantityConstraint(1,0,Inf));
     ump.addConstraint(QuantityConstraint(2,-Inf,0));
     ump.addConstraint(QuantityConstraint(3,-Inf,0));
     ump.addConstraint(QuantityConstraint(4,0,Inf));
%     for idx = 1:length(ump.instruments)
%         instrument=ump.instruments{idx};
%         ump.addConstraint(QuantityConstraint(idx,-instrument.bidSize,instrument.askSize));
%     end 
      [indifferencePrice,quantities]= ump.indifferencePrice(zcb,quantity, call, quantity*100);
%     plotPortfolio( ump.getInstruments(), quantities)
     [sup,sub] = ump.superHedgePrice( call,abs(quantity) );
end