classdef Future2 < Instrument1D
    %Future Represents holding a future. This class assumes a fixed
    % risk free rate so the bid and ask can be obtained by discounting
    
    properties
        %K is bid price for selling and ask price for buying
        %bid and ask are always zeros
        K
    end
    
    properties (Constant)
      contractSize=50;
      commission=0;
    end    
    methods
        function o = Future2(K,bid,ask,bidSize,askSize)
            % Construct the future given the quoted bid and ask as seen
            % on e.g. a Bloomberg terminal. 
            o = o@Instrument1D(bid,ask,bidSize,askSize);
            o.K=K;
        end
        
        % Compute the payoff of an instrument in the given scenarios.
        % each row of the scenarios matrix represnts a specific scenario
        % so this method should return a column vector
        function value = payoff( o, scenarios) 
            
            value = o.contractSize*(scenarios-o.K);

        end
        
        function K = getStrike(o)
            K = o.K;
        end
        
        function name = print( instrument)        
            name = sprintf('Future2 K=%d',instrument.K);
        end
        
        function delta = deltaAtInfinity(~)
            delta = 1;
        end        
        
        function wayPoints = getWaypoints(~)
            wayPoints = [];
        end
        
        
    end
    
end

