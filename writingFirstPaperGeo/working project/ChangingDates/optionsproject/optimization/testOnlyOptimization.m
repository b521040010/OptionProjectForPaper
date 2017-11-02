function utility=testOnlyOptimization()

    %delete(findall(0,'Type','figure'));
    initutilityoptimization();
    tic
    %riskAversion = 1*10^(-7);
    riskAversion = 0.00002;
    utilityFunction = ExponentialUtilityFunction( riskAversion );
    date = '20160408T145500';
    dayData = DayData( date );
  %  model = dayData.blackScholesModel();
     model= dayData.studentTModel();
    model.T=50/252;
    %BS model
    %For BS, we dont need to put log(S0) in there since we have a function
    %called logNormalParameters which will include log(S0) in mu later

   %  model.sigma=0.0713045/sqrt(model.T);
   %  model.mu=0.011272/model.T+0.5*model.sigma^2;
   %  model

    %Student-T model
    model.sigma=0.0553835;
    %model.sigma=0.4;
    model.mu=log(model.S0);
    model.nu=4.83548;
    %model.nu=20;
    arbitrage = ArbitrageFinder.findArbitrageForDate( date, false, true );
    assert(~arbitrage);
    % Create a utility maximization problem corresponding
    % to this problem
    ump = UtilityMaximizationProblem1D();
    ump.setModel( model );
    ump.setUtilityFunction(utilityFunction);
    % Set initial wealth
    zcb = dayData.findInstrument(0, DayData.cashType );
    
    for i=1:length(dayData.instruments)
        ump.addInstrument( dayData.instruments{i} );
    end    

%     K=2050;
%     try
%         call = dayData.findInstrument( K, DayData.callType );       
%         ump.removeInstrument( call );
%     catch
%         call=CallOption(K,1,1,1,1)
%     end
%     
%     try
%         put = dayData.findInstrument( K, DayData.putType );
%         ump.removeInstrument( put );
%     catch
%         put=PutOption(K,1,1,1,1)
%     end
%digital = DigitalCallOption(2050, 1,1,1000,1000);
%para=Parabola(2050, 1,1,1000,1000);
varianceSwap=VarianceSwap(2050, 1,1,1000,1000);
constant = ConstantForward(2050, 1,1,1000,1000);
    currentPort=Portfolio();
    %currentPort.add([100000],{varianceSwap});
    currentPort.add([100000+9.991684942319989e+03],{zcb});
    currentPort.add([-1],{constant});
    ump.setCurrentPosition(currentPort);
    
%      ump.addConstraint(QuantityConstraint(1,0,10^20));
%      ump.addConstraint(QuantityConstraint(2,-10^20,0));
%      ump.addConstraint(QuantityConstraint(3,-10^20,0));
%      ump.addConstraint(QuantityConstraint(4,0,10^20));
%    add constraints
% for idx = 1:length(ump.instruments)
%                 instrument=ump.instruments{idx};
%                  if isfinite(abs(instrument.bidSize))&&isfinite(abs(instrument.askSize))
%                     ump.addConstraint(QuantityConstraint(idx,-instrument.bidSize,instrument.askSize));
%                  end
%                  if isfinite(abs(instrument.bidSize))&&~isfinite(abs(instrument.askSize))
%                     ump.addConstraint(QuantityConstraint(idx,-instrument.bidSize,10^20));
%                  end
%                  if ~isfinite(abs(instrument.bidSize))&&isfinite(abs(instrument.askSize))
%                     ump.addConstraint(QuantityConstraint(idx,-10^20,instrument.askSize));
%                  end
%              end   
             for idx = 1:length(ump.instruments)
                    instrument=ump.instruments{idx};
                    ump.addConstraint(QuantityConstraint(idx,-instrument.bidSize,instrument.askSize));
             end      
    [utility, quantities] = ump.optimize();
   
    utility
    initialWealth=100000+9.991684942319989e+03;
    hold on
    ump.plotPortfolio( sprintf('Net profit of the portfolio invested on %s. Utility=%d', date,utility),quantities,initialWealth);
    plotPortfolio( ump.getInstruments(), quantities,initialWealth);
   % p=ump.returnPortfolio(quantities,initialWealth);
   % p.add([100000],{zcb})
%     
%     %Test if buying and selling quantities are in [-bidSizes,askSizes]
%     for idx=1:length(ump.instruments)
%         instrument=ump.instruments{idx};
%         assert(quantities(idx)<=instrument.askSize);
%         assert(quantities(idx)>=-instrument.bidSize);
%     end
    totalInvestment=0;
    for i = 1:length(ump.instruments)
        if quantities(i)>=0
            totalInvestment = totalInvestment+quantities(i)*ump.instruments{i}.getAsk();
        else
            totalInvestment = totalInvestment+quantities(i)*ump.instruments{i}.getBid();
        end
    end
    totalInvestment
%     temp=values(ump.currentPosition.map);
%     initialInvestment=temp{1}.quantity;
%     % Test that the initial invested money is less than our initial wealth
%     assert(totalInvestment<=initialInvestment);
end
