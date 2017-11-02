function [indifferencePrice,quantities]=testComputeIndifferencePriceWithoutLineSearch(riskAversion,K,quantity)

    delete(findall(0,'Type','figure'));
    initutilityoptimization();
    utilityFunction = ExponentialUtilityFunction( riskAversion );
    date = '20160408T145500';
    dayData = DayData( date );    
    arbitrage = ArbitrageFinder.findArbitrageForDate( date, false, true );
    assert(~arbitrage);

    %model = dayData.studentTModel();
    model = dayData.blackScholesModel();
    model.mu = 0.08;
    
    % Create a utility maximization problem corresponding
    % to this problem
    umpBarC = UtilityMaximizationProblem1D();
    umpBarCPlusC = UtilityMaximizationProblem1D();
    umpBarC.setModel( model );
    umpBarCPlusC.setModel( model );
    umpBarC.setUtilityFunction(utilityFunction);
    umpBarCPlusC.setUtilityFunction(utilityFunction);
    zcb = dayData.findInstrument(0, DayData.cashType );
    %my own adjustment from %ump.setCurrentPosition(zcb);
    currentPort=Portfolio();
    currentPort.add([100000],{zcb})
    umpBarC.setCurrentPosition(currentPort);    
    %umpBarCPlusC.setCurrentPosition(currentPort);  
    %------------------------
    currentPortBarC=Portfolio();
    currentPortBarC.add([100000],{zcb})
   

    for i=1:length(dayData.instruments)
        umpBarC.addInstrument( dayData.instruments{i} );
        umpBarCPlusC.addInstrument( dayData.instruments{i} );
    end
    
    umpBarC.addConstraint( BoundedLiabilityConstraint());
    umpBarCPlusC.addConstraint( BoundedLiabilityConstraint());

%      digitalCall=DigitalCallOption(K,1, 1.2 ,10,10)
    call = dayData.findInstrument( K, DayData.callType );
    put = dayData.findInstrument( K, DayData.putType );
    umpBarC.removeInstrument( call );
    umpBarC.removeInstrument( put );
    umpBarCPlusC.removeInstrument( call );
    umpBarCPlusC.removeInstrument( put );
     
    %add C to BarC
    currentPortBarC.add([quantity],{call})
    umpBarCPlusC.setCurrentPosition(currentPortBarC);  
    
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
% 
%      ump.addConstraint(QuantityConstraint(1,-Inf,Inf));
%     %ump.addConstraint(QuantityConstraint(2,-Inf,0));
%      ump.addConstraint(QuantityConstraint(2,-Inf,0));
%      ump.addConstraint(QuantityConstraint(3,0,Inf));


    for idx = 1:length(umpBarC.instruments)
        instrument=umpBarC.instruments{idx};
        umpBarC.addConstraint(QuantityConstraint(idx,-instrument.bidSize,instrument.askSize));
        umpBarCPlusC.addConstraint(QuantityConstraint(idx,-instrument.bidSize,instrument.askSize));
    end 
    
    [utilityBarC, ~] = umpBarC.optimize();
    [utilityBarCPlusC, ~] = umpBarCPlusC.optimize();
    
    utilityBarC=-(utilityBarC-1/riskAversion);
    utilityBarCPlusC=-(utilityBarCPlusC-1/riskAversion);
    indifferencePrice=(1/riskAversion)*log(utilityBarC/utilityBarCPlusC);
    quantities=[];
    
    
    
    
end