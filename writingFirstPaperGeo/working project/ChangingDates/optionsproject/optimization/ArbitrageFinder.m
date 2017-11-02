classdef ArbitrageFinder
    %ARBITRAGEFINDER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dayData
        ddd
    end
    
    methods
        
        function o = ArbitrageFinder( dayData )
%             for i =1:length(dayData.instruments)
%                 dayData.instruments{i}.bid=dayData.instruments{i}.bid/dayData.instruments{i}.contractSize;
%                 dayData.instruments{i}.bid=dayData.instruments{i}.ask/dayData.instruments{i}.contractSize;
%                 
%             end
            o.dayData = dayData;
            
            o.ddd = DoubledDayData(dayData);
        end
        
        function simpleInterest = getSimpleInterest(o)
            r = o.dayData.getInterestRate();
            T = o.dayData.getDaysToExpiry()/365.0;  
            simpleInterest = exp( r * T ) -1;
        end

        function [a, p] = findArbitrage( o, allowBonds, allowFutures )
            
            [a,p] = ArbitrageFinder2.findArbitrage( o.dayData.instruments, allowBonds, allowFutures );
            return;
        end        
        
    end
    
    methods (Static)

        function [arbitrage, portfolio] = findArbitrageForDate( date, allowBonds, allowFutures )
            dd = DayData( date );
            af = ArbitrageFinder( dd );
            [arbitrage, portfolio] = af.findArbitrage(allowBonds, allowFutures);
        end
        
        function findAllArbitrage(allowBonds, allowFutures)
            
            dayFiles = dir('../SPXFuturesAndOptions/');
            count = 0;
            arbitrageCount = 0;
            for i=1:length(dayFiles)
                file = dayFiles(i);
                name = file.name;
                if length(name)>4
                    date = name(1:end-4);
                    [foundArbitrage, portfolio] = ArbitrageFinder.findArbitrageForDate( date, allowBonds, allowFutures );
                    if foundArbitrage
                        disp( portfolio.print());
                    end
                    arbitrageCount = arbitrageCount + foundArbitrage;
                    fprintf('Date %s, found arbitrage in %f percent of files\n', date, arbitrageCount/count*100 );
                    count = count + 1;
                end
            end            
        end
                
        
    end
    
end

