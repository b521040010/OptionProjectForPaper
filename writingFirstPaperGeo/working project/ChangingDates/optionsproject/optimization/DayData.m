classdef DayData < matlab.mixin.Copyable
       
    
    properties
        spot;
        dateNum;                           
        instruments;
    end
    
    properties (Constant)
      cashType = 0;
      futureType = 1;
      callType = 2;
      putType = 3;
   end
    
    methods
        function [ dayData ] = DayData( dateString )

            % Reads market data for a particular day
            onlyDate = strcat(dateString(1:4),'-',dateString(5:6),'-',dateString(7:8))     
            dayData.dateNum = datenum(onlyDate,'yyyy-mm-dd');
            fileName = DayData.getFile( dateString );
            fprintf('Opening file %s\n', fileName );   
            matrix = csvread(fileName,1,1);
            dayData.spot=matrix(1,3);
            %dayData.spot=(matrix(1,3)+matrix(1,2))*0.5;
            % Convert to structured format
            strikes = matrix(:,1)';
            bids = matrix(:,2)';
            asks = matrix(:,3)';
            assert( sum(bids>asks)==0);
            isPuts = (matrix(:,4)==1)';
            isFuture = (matrix(:,5)==1)';
            bidSizes=matrix(:,6)';
            askSizes=matrix(:,7)';

            T = dayData.getT();
            r = dayData.getInterestRate();
            dayData.instruments = cell(0,1);
%             dayData.addInstrument( Bond( T, r, 1, 1,Inf,Inf ));
            %buy=paymoney
             
            %lending
            dayData.addInstrument( Bond( T, r, 1,1,0,Inf));
            %borrowing
            dayData.addInstrument( Bond( T, 0.03, 1, 1,Inf,0 ));
            %dayData.addInstrument( Bond( T, r, 1, 1,Inf,0 ));
            
            for i=2:length(bids)
                if isFuture(i)
                    dayData.addInstrument( Future2(bids(i),0,0,bidSizes(i),0));
                    dayData.addInstrument( Future2(asks(i),0,0,0,askSizes(i)));
                elseif isPuts(i)
                    dayData.addInstrument( PutOption(strikes(i),bids(i),asks(i),bidSizes(i),askSizes(i)));
                else
                    assert( strikes(i)>0 );                    
                    dayData.addInstrument( CallOption(strikes(i),bids(i),asks(i),bidSizes(i),askSizes(i)));
                end
            end            
        end        
        
        function price = askPrice(dayData, strike, instrumentType)
            ins = dayData.findInstrument(strike,instrumentType);
            price = ins.getAsk();
        end
        
        function price = bidPrice(dayData, strike, instrumentType)
            ins = dayData.findInstrument(strike,instrumentType);
            price = ins.getBid();
        end        
        
        function price = bidSize(dayData, strike, instrumentType)
            ins = dayData.findInstrument(strike,instrumentType);
            price = ins.getBidSize();
        end   
        
        function price = askSize(dayData, strike, instrumentType)
            ins = dayData.findInstrument(strike,instrumentType);
            price = ins.getAskSize();
        end     

        function price = contractSize(dayData, strike, instrumentType)
            ins = dayData.findInstrument(strike,instrumentType);
            price = ins.getContractSize();
        end             
        
        function price = commission(dayData, strike, instrumentType)
            ins = dayData.findInstrument(strike,instrumentType);
            price = ins.getCommission();
        end  
        
        function exists = containsInstrument(dayData, strike, instrumentType)
            [~, exists] = findInstrumentPrivate(dayData, strike, instrumentType);
        end        
        
        function ins = findInstrument(dayData, strike, instrumentType)
            [ins, exists] = findInstrumentPrivate(dayData, strike, instrumentType);
            if ~exists
                error('Instrument not found strike=%d, type=%d', strike, instrumentType);
            end
        end
        
        function [ins, exists] = findInstrumentPrivate(dayData, strike, instrumentType)
            exists = true;
            ins = [];
            n = length(dayData.instruments);
            for i=1:n
                instrument = dayData.instruments{i};
                if (instrumentType==DayData.cashType)
                    if isa(instrument,'Bond')
                        ins = instrument;
                        return;
                    end
                elseif (instrumentType==DayData.futureType)
                    if isa(instrument,'Future2') && instrument.getStrike()==strike
                        ins = instrument;
                        return;
                    end
                    if isa(instrument,'Future') && instrument.getStrike()==strike
                        ins = instrument;
                        return;
                    end                    
                elseif (instrumentType==DayData.callType)
                    if isa(instrument,'CallOption') && instrument.getStrike()==strike
                        ins = instrument;
                        return;
                    end
                elseif (instrumentType==DayData.putType)
                    if isa(instrument,'PutOption') && instrument.getStrike()==strike
                        ins = instrument;
                        return;
                    end                    
                end
            end
            exists = false;
        end
        
                    
        function interest = getInterestRate(dayData)
            interest = MarketData.getInstance().getInterestRate( dayData.dateNum );
        end
        
        function daysToExpiry = getDaysToExpiry(o)
            expiryDateNum = datenum(SPX.findNextExpiryDate(datevec(o.dateNum), 1));
            daysToExpiry = expiryDateNum - o.dateNum;
        end
        
        function T = getT(o)
            T = o.getDaysToExpiry()/365;
        end
        
        function m = blackScholesModel(dd)
            % Create a blackScholesModel calibrated to this day
            m = BlackScholesModel();
            m.T = dd.getT();
            m.S0 = dd.getSpot();
            
            nearStrike = round( m.S0/5 )*5;
            callPrice = dd.askPrice( nearStrike, DayData.callType );
            contractSize= dd.contractSize( nearStrike, DayData.callType );
            commission= dd.commission( nearStrike, DayData.callType );
            m.sigma = m.impliedVolatility( dd.getInterestRate(), nearStrike, false, (callPrice-commission)/contractSize);
        end
        
        function m=studentTModel(dd)
            m = StudentTModel();
            returns=m.getReturns;
            fitModel=m.fit( dd.getSpot(),dd.getT(), returns, ones(size(returns)));
            m.T=dd.getT();
            m.S0 = dd.getSpot();
            m.nu=fitModel.nu;
            m.sigma=sqrt(365)*fitModel.sigma;
            
            
        
        end
        
        function spot = getSpot(dd)
            spot=dd.spot;
        end
        
        function m = calibrateHistoric(dd, model, lambda )
            marketData = MarketData.getInstance();
            m = marketData.calibrateHistoric(model, ...
                                                      lambda, ...
                                                      dd.dateNum, ...
                                                      dd.spot,...
                                                      dd.getDaysToExpiry() );
        end
        
        function clearInstruments(o)
            o.instruments = cell(0,1);
        end
        
        function addInstrument( o, instrument ) 
            n = length(o.instruments);
            o.instruments{ n+1 } = instrument;
        end
        
    end
    
    methods (Static)
        
        function exists = dataExists( dateString )
            exists = exist( DayData.getFile( dateString ), 'file');
        end
        
        function fileName = getFile( dateString )            
            fileName = strcat( '../SPXFuturesAndOptions/',dateString,'.csv');
        end
        
        
    end
    
    

end

