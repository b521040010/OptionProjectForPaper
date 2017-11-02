classdef Portfolio < Instrument
    %PORTFOLIO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        map
    end
    
    methods
        function p = Portfolio()
            p@Instrument(NaN,NaN,Inf,Inf);
            p.map = containers.Map();
        end
                
        function add( p, quantities, instruments )
            assert( size( quantities,2)==1);
            assert( size( instruments,1)==1);
            assert( size( quantities, 1)==size(instruments,2));
            n = size(quantities,1);
            for i=1:n
                instrument = instruments{i};
                text = instrument.print();
                if (isKey(p.map,text))
                    entry = p.map(text);
                    entry.quantity = entry.quantity + quantities(i);
                    if (abs(entry.quantity)>=0)
                        p.map(text)=entry;
                    else
                        remove(p.map,text);
                    end
                else
                    entry.quantity = quantities(i);
                    entry.instrument = instruments{i};                    
                    if (abs(entry.quantity)>=0)
                        p.map(text)=entry;
                    end
                end
            end
        end
        
        function newPortfolio = updateMaturity( p, timePassedDays )
            newPortfolio = Portfolio();
            n = p.map.Count;
            ks = keys(p.map);
            for idx=1:n
                name = ks{idx};
                entry = p.map(name);
                q(1) = entry.quantity;
                ins{1} = entry.instrument.updateMaturity( timePassedDays );
                newPortfolio.add( q,ins);
            end
        end        
        
        function value = payoff( p, scenarios )
            n = p.map.Count;
            ks = keys(p.map);
            totalPayoff = zeros( size( scenarios ));
            for idx=1:n
                name = ks{idx};
                entry = p.map(name);
                iPayoff = entry.instrument.payoff( scenarios );
                totalPayoff = totalPayoff + iPayoff * entry.quantity;
            end
            value = totalPayoff;
        end
        
        function wayPoints = getWaypoints(p)
            n = p.map.Count;
            ks = keys(p.map);
            wayPoints = [];
            for idx=1:n
                name = ks{idx};
                entry = p.map(name);
                wayPoints = horzcat( wayPoints, entry.instrument.getWaypoints());
            end
            wayPoints = unique( wayPoints );
        end
        
        function d = deltaAtInfinity(p)
            % Compute the delta at infinity of the portfolio
            n = p.map.Count;
            ks = keys(p.map);
            d = 0;
            for idx=1:n
                name = ks{idx};
                entry = p.map(name);
                d = d + entry.quantity*entry.instrument.deltaAtInfinity();
            end
        end        
        
        % Print out the instrument returning a string
        function name = print( p ) 
            ks = keys(p.map);
            names = cell( p.map.Count+1, 1 );
            names{1}=sprintf('Portfolio:\n');
            for idx=1:p.map.Count
                name = ks{idx};
                entry = p.map(name);
                ins = entry.instrument;
                names{idx+1} = sprintf('%d * %s [bid=%d, ask=%d]\n', entry.quantity, name, ins.getBid(), ins.getAsk() );
            end
            name = sprintf('%s',names{:});
        end
        
        function markToMarket=computeMarkToMarket(p)
            sum=0;
            temp=p.map;
            temp=values(temp);
            for i=1:length(temp)
                instrument=temp{i}.instrument;
                mid=temp{i}.quantity*(instrument.bid+instrument.ask)/2;
                if instrument.contractSize==50
                    if temp{i}.quantity>=0
                        mid=(0.5*(instrument.bid+instrument.ask)-instrument.K);
                    else 
                        mid=(0.5*(instrument.bid+instrument.ask)-instrument.K);
                    end
                end
                sum=sum+mid;
            end
            markToMarket=sum;
        end
    end
    
end

