classdef OldPortfolio 
    properties
        % Vector of strikes
        strikes
        % Associated vector of types
        instrumentTypes
        % Associated vector of quantities of that option
        quantities
    end
    methods        
        
        function p = OldPortfolio()
            % Create a simple example portfolio
            p.quantities = [-0.5, 0.8, -0.3];
            p.strikes = [0, 1.0,1.1];
            p.instrumentTypes = [ DayData.cashType, DayData.callType , DayData.putType];            
        end
        
        function r = payoff( p, simpleInterest, askPrice )
            % According to SPX definitions, the options payoff at the
            % "opening sales price" - we'll take this to be the ask price
            nInstruments = length(p.strikes);
            nX = length(askPrice);
            askMat = repmat(askPrice',[1,nInstruments]);
            Lstrikes = repmat(p.strikes,[nX,1]);
            
            isCalls = p.instrumentTypes == DayData.callType;
            isPuts = p.instrumentTypes == DayData.putType;
            isCash = p.instrumentTypes == DayData.cashType;            
            isStock = p.instrumentTypes == DayData.futureType;
            
            LisCalls = repmat(isCalls,[nX,1]);
            LisPuts = repmat(isPuts,[nX,1]);
            LisStock = repmat(isStock,[nX,1]);
            Lquantities = repmat(p.quantities,[nX,1]);
            
            cash = sum(p.quantities(isCash));
            payoffs = positivePart(askMat - Lstrikes).* LisCalls ...
                    + positivePart(Lstrikes-askMat).* LisPuts ...
                    + (askMat-Lstrikes) .* LisStock;
            r = sum( payoffs .* Lquantities, 2 ) + (1.0 + simpleInterest)*cash;                        
        end
        
        function plotPayoff( p, title, simpleInterest )
            figure('Name',title,'NumberTitle','off');
            to = 1.1*max(p.strikes);
            if to==0
                to = 100;
            end
            x = linspace(0,to,1000);
            y = p.payoff( simpleInterest, x, x );
            plot(x,y);
        end
        
        function p = simplify( p )
            % Eliminate unnecessary duplicates from a portfolio
            
            newStrikes = [];
            newQuantities = [];
            newTypes = [];
            
            for type=[DayData.putType DayData.callType DayData.futureType DayData.cashType]
                quantityMap = containers.Map('KeyType','double','ValueType','double');
                matchingStrikes = p.strikes( p.instrumentTypes==type );
                matchingQuantities = p.quantities( p.instrumentTypes==type );
                for i=1:length(matchingStrikes)
                    key = matchingStrikes(i);
                    if isKey(quantityMap,key)
                        quantityMap(key) = quantityMap(key) + matchingQuantities(i);
                    else
                        quantityMap(key) =  matchingQuantities(i);
                    end
                end
                newStrikes = [cell2mat(keys( quantityMap )) newStrikes ];
                newQuantities = [cell2mat(values( quantityMap )) newQuantities ];
                newTypes = [type*ones(1, quantityMap.Count) newTypes ];
            end
            
            p.strikes = newStrikes;
            p.quantities = newQuantities;
            p.instrumentTypes = newTypes;
        end                
        
        function total = cost( p, dayData )
            total= marketPrice(p, dayData,1);
        end                

        function total = value( p, dayData )
            total= marketPrice(p, dayData,-1);
        end                
        
        function total = marketPrice( p, dayData, sgn )
            % How much would it cost to purchase the portfolio?
            total = 0;
            nInstruments = length( p.strikes )
            
            for i=1:nInstruments
                askFilter = ((sgn*p.quantities(i))>0)
                bidFilter = ((sgn*p.quantities(i))<0)
                if max(askFilter)>0
                    %%%%%%%%%%%%%%%%%%%%%%%%%%
                    if sgn==1
                        total = total ...
                            + dayData.askPrice(p.strikes(i),p.instrumentTypes(i)) .* p.quantities(i) .* askFilter;
                    else
                        dayData.instruments{i}
                        dayData.instruments{i}.commission
                        total = total ...
                            + dayData.askPrice(p.strikes(i),p.instrumentTypes(i)) .* p.quantities(i) .* askFilter-2*dayData.instruments{i}.commission.* p.quantities(i);
                    end
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         total = total ...
%                             + dayData.askPrice(p.strikes(i),p.instrumentTypes(i)) .* p.quantities(i) .* askFilter;                    
                end                    
                if max(bidFilter)>0
                    %%%%%%%%%%%%%%%%%%%%%%%%%%
                     if sgn==1
                         total = total ...
                        + dayData.bidPrice(p.strikes(i),p.instrumentTypes(i)) .* p.quantities(i) .* (bidFilter);
                     else
                         dayData.instruments{i}
                         dayData.instruments{i}.commission
                         total = total ...
                        + dayData.bidPrice(p.strikes(i),p.instrumentTypes(i)) .* p.quantities(i) .* (bidFilter)+2*dayData.instruments{i}.commission.* p.quantities(i);
                     end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     total = total ...
%                         + dayData.bidPrice(p.strikes(i),p.instrumentTypes(i)) .* p.quantities(i) .* (bidFilter);
                end                    
            end
        end            
        
        function q = quantity( p, strike, instrumentType )
            % Returns the quantity of the given instrument
            q = p.quantities( p.strikes==strike & p.instrumentTypes==instrumentType);
        end

        function disp( p )
            fprintf('Portfolio\n');
            nInstruments = length(p.strikes);
            for i=1:nInstruments
                if p.quantities(i)~=0
                    fprintf('Strike %f, type %f, quantity %f\n', p.strikes(i),p.instrumentTypes(i),p.quantities(i));
                end        
            end
        end
        
    end
end

