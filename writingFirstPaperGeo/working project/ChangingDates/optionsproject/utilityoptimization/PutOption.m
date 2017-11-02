classdef PutOption < Instrument1D
    %PUTOPTION Represents a put option
    
    properties
        K
    end
    
    properties (Constant)
      contractSize=100;
      commission=0;
      %commission=0;
    end
    
    methods
        
        function o = PutOption(K, bid,ask,bidSize,askSize)
            o = o@Instrument1D(bid,ask,bidSize,askSize);
            %o.bid=o.bid*o.contractSize;
            %o.ask=o.ask*o.contractSize;
            o.bid=o.bid*o.contractSize-o.commission;
            o.ask=o.ask*o.contractSize+o.commission;              
            o.K = K;
        end        
        
        
        % Compute the payoff of an instrument in the given scenarios.
        % each row of the scenarios matrix represnts a specific scenario
        % so this method should return a column vector
        function value = payoff( instrument, scenarios ) 
            diff = instrument.K-scenarios;
            value = instrument.contractSize*(diff > 0) .* diff;
        end
        
        function name = print( instrument)        
            name = sprintf('Put K=%d',instrument.K);
        end
        
        function delta = deltaAtInfinity(~)
            delta = 0;
        end        
        
        function wayPoints = getWaypoints(o)
            wayPoints = o.K;
        end
        
        function K = getStrike(o)
            K = o.K;
        end     

        function K = getCommission(o)
            K = o.commission;
        end 
        
        function K = getContractSize(o)
            K = o.contractSize;
        end             
    end
    
end

