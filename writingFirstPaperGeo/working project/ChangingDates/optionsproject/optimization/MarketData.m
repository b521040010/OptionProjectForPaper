classdef MarketData
    
    properties
        spxIndex;
        spxSettlementIndex;
        us0001mIndex;
    end
    
    methods(Static)
        function ret = getInstance( reset )
            % Access to singleton instance
            persistent marketData;
            persistent marketDataInitialized;            
            
            if nargin<1
                reset = false;
            end
                        
            if isempty(marketDataInitialized) || reset
                marketData = MarketData();
                marketDataInitialized = true;
            end
            ret = marketData;
        end
    end
    
    methods (Access=private)
        function marketData = MarketData()
            disp('Initializing market data');
            % Constructor 
            spxI.values = xlsread('../SPX Index.xlsx', 1, 'B3:B2519');
            [~,dateStrings] =  xlsread('../SPX Index.xlsx', 1, 'A3:A2519');
            spxI.dates = datenum( dateStrings, 'dd/mm/yyyy');

            spxSetI.values = xlsread('../SPXSettlement.xlsx', 1, 'B3:B1009');
            [~,dateStrings] =  xlsread('../SPXSettlement.xlsx', 1, 'A3:A1009');
            spxSetI.dates = datenum( dateStrings, 'dd/mm/yyyy');
            
            us1MI.values = xlsread('../US0001MIndex.xlsx', 1, 'D3:D1034');
            [~,dateStrings] =  xlsread('../US0001MIndex.xlsx', 1, 'A3:A1034');
            us1MI.dates = datenum( dateStrings, 'dd/mm/yyyy');

            marketData.spxIndex = spxI;
            marketData.spxSettlementIndex = spxSetI;
            marketData.us0001mIndex = us1MI;     
            disp('Market data initialized');
        end
    end
    
    methods
        
        function model = calibrateHistoric( md, model, lambda, dateNum, S0, daysToExpiry )
            % Calibrate the model using historic returns data
            % Lambda is the exponential weighting parameter
            numberOfDays = daysToExpiry;
            assert( numberOfDays>=1 );
            idx = find( md.spxIndex.dates <= dateNum, 1, 'last');
            returns = []; 
            count = 1;
            keepLooping = ~isempty(idx);
            while keepLooping
                date = md.spxIndex.dates(idx);
                value = md.spxIndex.values(idx);
                prevIdx = find( md.spxIndex.dates < date-numberOfDays, 1, 'last');
                if (idx==prevIdx+1)
                    prevIdx = prevIdx-1;
                end
                if (~isempty(prevIdx))
              
      prevDate = md.spxIndex.dates(prevIdx+1);
                    prevValue = md.spxIndex.values(prevIdx+1);
                    actualNDays = date-prevDate;
                    actualReturn = (value-prevValue)/prevValue;
                    scaledReturn = actualReturn * sqrt( numberOfDays/actualNDays);
                    returns(count) = scaledReturn;
                    count = count + 1;
                    idx = prevIdx+1;
                else
                    keepLooping = 0;
                end
            end
            
            weights = lambda .^ (1:length(returns)) * 1e6;
            weights = round(weights);
            returns = returns(weights>0);
            weights = weights(weights>0);            
            
            model = model.fit( S0, numberOfDays/365, returns', weights');
        end
        
        % Returns the annualized interest rate
        function ir = getInterestRate( marketData, dateNum )
            us1MI = marketData.us0001mIndex;
            idx = find( us1MI.dates >= dateNum, 1, 'first');
            if us1MI.dates(idx)==dateNum
                ir = us1MI.values(idx);
            else
                assert( idx>1 );
                d1 = abs(us1MI.dates(idx-1)-dateNum);
                d2 = abs(us1MI.dates(idx)-dateNum);
                weight1 = d2/(d1+d2);
                weight2 = d1/(d1+d2);
                ir = us1MI.values(idx-1)*weight1 + us1MI.values(idx)*weight2;
            end     
            ir = log(1+ir/100);
        end
        
        % Returns the settlement value on the given date
        function value = getSPXSettlementValue( marketData, dateNum )
            spxSI = marketData.spxSettlementIndex;
            idx = find( spxSI.dates == dateNum, 1, 'first');
            if spxSI.dates(idx)==dateNum
                value = spxSI.values(idx);
            else
                error( 'Settlement value not available for %d', dateNum );
            end     
        end
        
        % Returns the SPX index value on the given date
        function value = getSPXIndexValue( marketData, dateNum )
            spxI = marketData.spxIndex;
            idx = find( spxI.dates == dateNum, 1, 'first');
            if spxI.dates(idx)==dateNum
                value = spxI.values(idx);
            else
                error( 'SPX value not available for %d', dateNum );
            end     
        end        
        
    end
    
end

