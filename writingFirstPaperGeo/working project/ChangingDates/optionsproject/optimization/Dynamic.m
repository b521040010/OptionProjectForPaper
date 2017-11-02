
classdef Dynamic < matlab.mixin.Copyable
    
    properties
        startDate
        endDate
        rebalancingInterval
        initialWealth
        histDiffPort
        histPort
        histQuantities
        histUtility
        currentState
        dateString
        histSpot
        quantities
        qp
        s
        vol
    end
    
    methods
        
        function o=Dynamic(startDate,endDate,rebalancingInterval,initialWealth)
            o.startDate=startDate;
            o.endDate=endDate;
            o.rebalancingInterval=rebalancingInterval;
            o.initialWealth=initialWealth;
            
                dayData = DayData( startDate );
                zcb = dayData.findInstrument(0, DayData.cashType );
                currentPort=Portfolio();
                currentPort.add([initialWealth],{zcb});
            o.histPort=struct('initialWealthhhh',currentPort);
            o.histDiffPort=struct('initialWealthhhh',currentPort);
            o.histUtility=struct('initialWealthhhh',0);
            o.histSpot=struct('initialWealthhhh',dayData.spot);
            o.vol=struct('initialWealthhhh',0);
            o.currentState=2;
        end
        
        function o=run(o)
            tic
            o=o.prepareTheDates;
            
            for i=2 : size(o.dateString,1)
                o=o.reoptimize;
                o.currentState=o.currentState+1;
             
               
            end 
            
            
        end
        
        function o=prepareTheDates(o)
            %If the investment horizon is 20 days, and rebalancingInterval
            %is 6 days, we will invest only 3 times. The last investment
            %horizon is 8 days.

            onlyDateStartDate = strcat(o.startDate(2:5),'-',o.startDate(6:7),'-',o.startDate(8:9));
            onlyDateEndDate = strcat(o.endDate(2:5),'-',o.endDate(6:7),'-',o.endDate(8:9));
            numOnlyDateStartDate=datenum(onlyDateStartDate);
            numOnlyDateEndDate=datenum(onlyDateEndDate);
            numberOfRebalancing=floor((numOnlyDateEndDate -numOnlyDateStartDate )/o.rebalancingInterval);
            temp=o.rebalancingInterval*[0 ones(1,numberOfRebalancing-1)];
            tempDate=numOnlyDateStartDate+cumsum(temp);
            dateTimeZero=datestr(tempDate,'yyyymmddTHHMMSS');
%             o.dateString=strcat('D',dateTimeZero(:,1:8),o.startDate(10:ek.ud));
%             o.dateString=vertcat('initialWealthhhh',o.dateString);
            
            
            j=2;
            o.dateString='initialWealthhhh';
            dateWithHoliday=strcat('D',dateTimeZero(:,1:8),o.startDate(10:end));
            for i =1:size(dateWithHoliday,1)
                try 
                    xlsread(strcat( '../SPXFuturesAndOptions/',dateWithHoliday(i,:),'.csv'));
                    o.dateString=vertcat(o.dateString,dateWithHoliday(i,:));
                    j=j+1;
                catch
                end
            end
            
            
            
        end
        
        function o=reoptimize(o)
            riskAversion = 0.00002;
            utilityFunction = ExponentialUtilityFunction( riskAversion );
            date = o.dateString(o.currentState,:);
            previousDate=o.dateString(o.currentState-1,:);
            dayData = DayData( date );
           model = dayData.blackScholesModel();
            
            model.mu = 0.5*model.sigma^2;
  %           model.mu=0;
%             model.sigma=0.05;
%            model
%                model= dayData.studentTModel();
%              model.mu=log(model.S0);
%              model.sigma=0.0553835;
%              model.mu=0.0173861+log(model.S0);
%              model.nu=4.83548;

            model
%             arbitrage = ArbitrageFinder.findArbitrageForDate( date, false, true );
%             assert(~arbitrage);    
            ump = UtilityMaximizationProblem1D();
            ump.setModel( model );
            ump.setUtilityFunction(utilityFunction);
            %currentPort=(o.histPort.(previousDate));
            ump.setCurrentPosition(o.histPort.(previousDate));
%            ump.setCurrentPosition(currentPort);

            for i=1:length(dayData.instruments)
                ump.addInstrument( dayData.instruments{i} );
            end    
            for idx = 1:length(ump.instruments)
                instrument=ump.instruments{idx};
                 if isfinite(abs(instrument.bidSize))&&isfinite(abs(instrument.askSize))
%                      if abs(instrument.bidSize)>50
%                          instrument.bidSize=50;
%                      end
%                      if abs(instrument.askSize)>50
%                          instrument.askSize=50;
%                      end
                ump.addConstraint(QuantityConstraint(idx,-instrument.bidSize,instrument.askSize));
                 end
%                 if isfinite(abs(instrument.bidSize))&&isfinite(abs(instrument.askSize))
%                 idx;
%                 ump.addConstraint(QuantityConstraint(idx,-50,51));
%                 end
             end    
            %ump.addConstraint( BoundedLiabilityConstraint());
            
%             prices = model.simulatePricePaths(1000000,1);
%             scenarios = prices(:,end);
%             o.s = sort(scenarios);
            
            [utility, quantities,qp] = ump.optimize();
            utility
            o.quantities=quantities;
            o.qp=qp;
            currentPort=(o.histPort.(previousDate));
            temp=values(currentPort.map);
            ii=cell(1,length(temp));
                for i =1:length(temp)
                    tempp=temp{i};
                    qq(i)=tempp.quantity;
                    ii{i}=tempp.instrument;
                end
            port=Portfolio();
            port.add(quantities,ump.instruments);
            port.add(qq',ii);
            
            
            diffPort=Portfolio();
            diffPort.add(quantities,ump.instruments);
            
                totalInvestment=0;
    for i = 1:length(ump.instruments)
        if quantities(i)>=0
%             i
%             ask=ump.instruments{i}.getAsk()
            totalInvestment = totalInvestment+quantities(i)*ump.instruments{i}.getAsk();
        else
%             i
%             bid=ump.instruments{i}.getBid()
            totalInvestment = totalInvestment+quantities(i)*ump.instruments{i}.getBid();
        end
    end
    totalInvestment
    
    newUtility=utility+(1-utility*0.00002)*(1-exp(0.00002*totalInvestment))/0.00002
     assert(totalInvestment<=5)
     assert(totalInvestment>=-5)
            
            o.histPort=setfield(o.histPort,date,port);
            o.histDiffPort=setfield(o.histDiffPort,date,diffPort);
            o.histUtility=setfield(o.histUtility,date,utility);
            o.histSpot=setfield(o.histSpot,date,dayData.spot);
            o.vol=setfield(o.vol,date,model.sigma);

            
                %Test if buying and selling quantities are in [-bidSizes,askSizes]
%     for idx=1:length(ump.instruments)
%         instrument=ump.instruments{idx};
%         assert(quantities(idx)<=instrument.askSize);
%         assert(quantities(idx)>=-instrument.bidSize);
%     end

            
            
            
        end
        
        function plotHistPortfolio(o,date)
            port=o.histPort.(date);
            temp=values(port.map);
            instruments=cell(length(temp),1);
            quantities=zeros(length(temp),1);
            for i=1:length(temp)
                instruments{i}=temp{i}.instrument;
                quantities(i)=temp{i}.quantity;
            end
            plotPortfolioWithoutFutures(date,instruments,quantities);
        end
        
        function plotHistDiffPortfolio(o,date)
            port=o.histDiffPort.(date);
            temp=values(port.map);
            instruments=cell(length(temp),1);
            quantities=zeros(length(temp),1);
            for i=1:length(temp)
                instruments{i}=temp{i}.instrument;
                quantities(i)=temp{i}.quantity;
            end
            plotPortfolioWithoutFutures(date,instruments,quantities);
        end
        
        function [adjustedMark adjustedSpot]=getMarkToMarket(o)
            date=o.dateString;
            mark=zeros(1,length(date));
            port= o.histPort.(date(1,:));
            mark(1)=port.computeMarkToMarket;
            spot(1)= o.histSpot.(date(1,:));
            for i=2:length(date)
                port=o.histPort.(date(i,:));
                mark(i)=port.computeMarkToMarket;
                spot(i)= o.histSpot.(date(i,:));
            end
            adjustedMark=mark(2:end);
            adjustedSpot=spot(2:end)

            
        end
        
        function plotMarkToMarket(o)
            date=o.dateString;
            mark=zeros(1,size(date,1));
            port= o.histPort.(date(1,:));
            mark(1)=port.computeMarkToMarket;
            spot(1)= o.histSpot.(date(1,:));
            for i=2:size(date,1)
                port=o.histPort.(date(i,:));
                mark(i)=port.computeMarkToMarket;
                spot(i)= o.histSpot.(date(i,:));
            end
            subplot(2,1,1);
            plot((1:1:length(mark)-1),mark(2:end));
            subplot(2,1,2);
            plot((1:1:length(mark)-1),spot(2:end));
            
        end

 function plotLogMarkToMarket(o)
            date=o.dateString;
            mark=zeros(1,length(date));
            port= o.histPort.(date(1,:));
            mark(1)=port.computeMarkToMarket;
            spot(1)= o.histSpot.(date(1,:));
            for i=2:length(date)
                port=o.histPort.(date(i,:));
                mark(i)=port.computeMarkToMarket;
                spot(i)= o.histSpot.(date(i,:));
            end
            subplot(2,1,1);
            plot((1:1:length(mark)-1),log(mark(2:end)));
            subplot(2,1,2);
            plot((1:1:length(mark)-1),spot(2:end));
            
        end        
        
        function plotMarkToMarketPercent(o)
            date=o.dateString;
            mark=zeros(1,length(date));
            port= o.histPort.(date(1,:));
            mark(1)=port.computeMarkToMarket;
            spot(1)= o.histSpot.(date(1,:));
            for i=2:length(date)
                port=o.histPort.(date(i,:));
                mark(i)=port.computeMarkToMarket;
                spot(i)= o.histSpot.(date(i,:));
            end
            subplot(2,1,1);
            plot((1:1:length(mark)-1),100*(mark(2:end)-o.initialWealth)/(o.initialWealth));
            subplot(2,1,2);
            plot((1:1:length(mark)-1),100*(spot(2:end)-spot(1))/spot(1));
            
        end
        
        function plotMarkToMarketDevelopePercent(o)
            date=o.dateString;
            mark=zeros(1,length(date));
            port= o.histPort.(date(1,:));
            mark(1)=port.computeMarkToMarket;
            spot(1)= o.histSpot.(date(1,:));
            for i=2:length(date)
                port=o.histPort.(date(i,:));
                mark(i)=port.computeMarkToMarket;
                spot(i)= o.histSpot.(date(i,:));
            end
            subplot(2,1,1);
            bar((1:1:length(mark)-1),100*(mark(2:end)-mark(1:end-1))./(mark(1:end-1)));
            subplot(2,1,2);
            bar((1:1:length(mark)-1),100*(spot(2:end)-spot(1:end-1))./spot(1:end-1));
            
        end
        
        function plotAll(o)
            o.plotMarkToMarket
            figure
            o.plotLogMarkToMarket
            figure
            o.plotMarkToMarketPercent
            figure
            o.plotMarkToMarketDevelopePercent
        
        end
         function plotLogLog(o)
            date=o.dateString;
            mark=zeros(1,length(date));
            port= o.histPort.(date(1,:));
            mark(1)=port.computeMarkToMarket;
            spot(1)= o.histSpot.(date(1,:));
            for i=2:length(date)
                port=o.histPort.(date(i,:));
                mark(i)=port.computeMarkToMarket;
                spot(i)= o.histSpot.(date(i,:));
            end
            plot((1:1:length(mark)-1),log(mark(2:end)));
            hold on
            plot((1:1:length(mark)-1),log(mark(2))-log(spot(2))+log(spot(2:end) ));
            
         end    
         function plotLogLogScaled(o)
            date=o.dateString;
            mark=zeros(1,length(date));
            port= o.histPort.(date(1,:));
            mark(1)=port.computeMarkToMarket;
            spot(1)= o.histSpot.(date(1,:));
            for i=2:length(date)
                port=o.histPort.(date(i,:));
                mark(i)=port.computeMarkToMarket;
                spot(i)= o.histSpot.(date(i,:));
            end

            plot((1:1:length(mark)-1),log(mark(2:end)/mark(2)));
            hold on
            plot((1:1:length(mark)-1),log(spot(2:end)/spot(2) ));
            
        end    
        
    end
end