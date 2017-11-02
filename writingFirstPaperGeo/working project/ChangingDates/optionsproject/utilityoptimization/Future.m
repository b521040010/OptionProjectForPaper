classdef Future < Instrument1D
    %Future Represents holding a future. This class assumes a fixed
    % risk free rate so the bid and ask can be obtained by discounting
    
    properties
    end
    
    properties (Constant)
        contractSize=1;
      commission=0;
    end
    
    methods
        function o = Future(bid,ask,bidSize,askSize)
            % Construct the future given the quoted bid and ask as seen
            % on e.g. a Bloomberg terminal. 
            o = o@Instrument1D(bid,ask,bidSize,askSize);
        end
        
        % Compute the payoff of an instrument in the given scenarios.
        % each row of the scenarios matrix represnts a specific scenario
        % so this method should return a column vector
        function value = payoff( o, scenarios ) 
            
            value = scenarios;
        end
        
        
        function name = print( instrument)        
            name = 'Future';
        end
        
        function delta = deltaAtInfinity(~)
            delta = 1;
        end        
        
        function wayPoints = getWaypoints(~)
            wayPoints = [];
        end
        
        
    end
    
end

