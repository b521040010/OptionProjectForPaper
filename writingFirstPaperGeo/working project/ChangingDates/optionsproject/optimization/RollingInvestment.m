classdef RollingInvestment < matlab.mixin.Copyable
    %RollingInvestment
    %   This simulates an investment strategy over a given time period
    
    properties
        startDate;
        endDate;
        
        currentDate;
        currentDayData;
        currentPortfolio;   
        currentPrice;
        
        currentInterestRate;
        
        days;
        superHedgeHistory;
        subHedgeHistory;
        priceHistory;
        spx;
    end
    
    methods
        
        function o = RollingInvestment( startDate, endDate, initialCapital )
            o.startDate = startDate;
            o.endDate = endDate;
            o.currentDate = startDate;
            o.currentPrice = initialCapital;
            o.currentDayData = DayData( o.currentDate );
            
            o.initBondPortfolio( initialCapital );
        end       
        
        function initBondPortfolio( o, value ) 
            o.currentInterestRate = o.currentDayData.getInterestRate();
            
            o.currentPortfolio = Portfolio();
            % Add an instrument corresponding to cash under the mattress
            q(1) = value;
            i{1} = Bond(1,o.currentInterestRate,1,1);
            o.currentPortfolio.add( q, i);            
        end
        
        function run( o ) 
            dateNumStart = datenum(o.startDate,'yyyy-mm-dd');
            dateNumEnd = datenum(o.endDate,'yyyy-mm-dd');
            for dateNum = dateNumStart:dateNumEnd
                dateString = datestr(dateNum,'yyyy-mm-dd','local');
                if (DayData.dataExists(dateString))
                    o.processDate( dateString );
                end
            end
        end

        function processDate( o, dateString )
            prevDate = o.currentDate;
            o.currentDate = dateString;
            timePassed = datenum( dateString ) - datenum( prevDate );
            o.currentPortfolio = updateMaturity( o.currentPortfolio, timePassed );
            
            o.currentDayData = DayData( dateString );
            daysToExpiry = o.currentDayData.getDaysToExpiry();            
            if daysToExpiry==1
                expiryDate = o.currentDayData.dateNum;
                settlementValue = MarketData.getInstance().getSPXSettlementValue( expiryDate );
                portfolioValue = o.currentPortfolio.payoff( settlementValue );
                o.recordValue( portfolioValue, portfolioValue, portfolioValue );
                o.initBondPortfolio( portfolioValue );
            else
                if o.futurePriceFound()
                    %arbitrage = ArbitrageFinder.findArbitrageForDate( o.currentDate, true, true );
                    %if (~arbitrage)
                        umpParams = o.createUmpParams();
                        ump = o.createUtilityMaximizationProblem(umpParams);
                        try 
                            [~, quantities] = ump.optimize();
                            instruments = ump.getInstruments();                        
                            newPortfolio = o.currentPortfolio;
                            newPortfolio.add( quantities, instruments );
                            o.currentPortfolio = newPortfolio;
                        catch ex1
                            disp('Optimization failure');
                        end
                        % price the portfolio
                        T = o.currentDayData.getDaysToExpiry()/365.0;
                        zcb = Bond( T, o.currentDayData.getInterestRate(), 1, 1 );
                        
                        umpParams = o.createUmpParams();
                        umpParams.nSdsForOptions = 2;
                        ump = o.createUtilityMaximizationProblem( umpParams);
                        ump.setCurrentPosition( Portfolio() );
                        ump.addInstrument( zcb );
                        try 
                            [sup, sub] = ump.superHedgePrice( o.currentPortfolio );                    
                        catch ex2
                            disp('Failed to calculate sup and sub hedge prices for at a time point');
                            sup = NaN;
                            sub = NaN;
                        end
                        try 
                            indifferencePrice = ump.indifferencePriceExponentialUtility(zcb, 1, o.currentPortfolio, o.currentPrice );
                        catch ex3
                            disp('Failed to calculate indifference price for a time point time point');
                            indifferencePrice = NaN;
                        end
                        o.recordValue( sub, sup, indifferencePrice );
                    %end
                end
            end
            o.plot();
        end
        
        function f = futurePriceFound( o )
            % Is there a future price for the given day?
            f = o.currentDayData.containsInstrument( 0, DayData.futureType );
        end
        
        function umpParams = createUmpParams(~)
            % Create parameters used to configure the utility maximization
            % problem
            umpParams.riskAversion=1.0;
            umpParams.ewmaLambda = 0.98;
            umpParams.nSdsForOptions = 10; % Only trade in options that are reasonably in the money
            umpParams.quantityConstraint = 10000;
            umpParams.model = StudentTModel();
        end
        
        
        function ump = createUtilityMaximizationProblem(o, umpParams)
            % Create utility maximization problem corresponding to the
            % parameters
                        
            utilityFunction = ExponentialUtilityFunction( umpParams.riskAversion );
            model = o.currentDayData.calibrateHistoric( umpParams.model, umpParams.ewmaLambda);
            expiryDays = o.currentDayData.getDaysToExpiry();
            ump = UtilityMaximizationProblem1D();
            ump.setModel( model );
            ump.setUtilityFunction( utilityFunction );

            ump.setCurrentPosition( o.currentPortfolio );

            for idx=1:length(o.currentDayData.instruments)
                ump.addInstrument( o.currentDayData.instruments{idx} );
            end  
            
            
            spot = o.currentDayData.getSpot();
            o.removeOutOfTheMoney( ump, spot, umpParams.model.sigma, expiryDays, umpParams.nSdsForOptions ); 
            % Add constraint on quantity
            instruments = ump.getInstruments();
            for idx=1:length(instruments)
                price = instruments{idx}.getAsk();
                if isnan(price)
                    price = instruments{idx}.getBid();
                end                
                qc = QuantityConstraint(idx,-umpParams.quantityConstraint/price,umpParams.quantityConstraint/price);
                ump.addConstraint( qc );
            end  

            
        end
        
        function removeOutOfTheMoney( ~, ump, spot, sigma, expiryDays, nSds) 
            min = spot - nSds*sigma*sqrt( expiryDays/365.0 )*spot;
            max = spot + nSds*sigma*sqrt( expiryDays/365.0 )*spot;
            function accept = filter( instrument )
                if (isa(instrument,'CallOption') || isa(instrument,'PutOption'))
                    accept = instrument.K>=min && instrument.K<=max;
                elseif (isa(instrument,'Bond'))
                    accept = 0;
                else
                    accept = 0;
                end
            end
            ump.filterInstruments( @filter );
        end
                
        % Record the current value of the portfolio
        function recordValue( o, sub, sup, currentPrice)                                                                                   
            fprintf('Price in [%d,%d]\n', sub, sup);
            n = length( o.superHedgeHistory );
            o.currentPrice = currentPrice;
            o.subHedgeHistory(n+1)= sub;
            o.superHedgeHistory(n+1)= sup;
            o.priceHistory(n+1)=currentPrice;
            md = MarketData.getInstance();
            o.spx(n+1) = md.getSPXIndexValue( datenum(o.currentDate) );
            o.days(n+1) = datenum(o.currentDate ) - datenum( o.startDate );                        
        end
        
        function plot(ri)
            md = MarketData.getInstance();
            spxStart = md.getSPXIndexValue(datenum(ri.startDate) );
            %plot( ri.days, ri.subHedgeHistory, '-r');
            plot( ri.days, ri.priceHistory, '-g');
            hold on;
            %plot( ri.days, ri.superHedgeHistory, '-b');
            plot( ri.days, ri.spx/spxStart, '-m' );
            %legend('Sub hedge value', 'Indifference Price', 'Super hedge value', 'SPX Index');
            legend('Indifference Price', 'SPX Index');
            hold off;            
        end
                               
    end
    
    
    
end



