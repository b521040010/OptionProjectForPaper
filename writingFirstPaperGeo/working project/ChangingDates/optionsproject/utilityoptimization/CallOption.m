classdef CallOption < Instrument1D
    %CALLOPTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        K
        
    end
    
    properties (Constant)
      contractSize=100;
      commission=0;
% commission is qouted in terms of percentage
    %   commission=0;
    end
    
    methods
        
        function o = CallOption(K, bid,ask,bidSize,askSize)
            o = o@Instrument1D(bid,ask,bidSize,askSize);
            %o.bid=o.bid*o.contractSize;
            %o.ask=o.ask*o.contractSize;
            o.bid=o.bid*o.contractSize-o.commission;
            o.ask=o.ask*o.contractSize+o.commission;            
            o.K = K;
        end        
        
        % Compute the payoff of an instrument in the given sceonarios.
        % each row of the scenarios matrix represnts a specific scenario
        % so this method should return a column vector
        function value = payoff( instrument, scenarios ) 
            diff = scenarios-instrument.K;
            value = instrument.contractSize*(diff > 0) .* diff;
        end
        
        function delta = deltaAtInfinity(~)
            delta = 1;
        end        
        
        function wayPoints = getWaypoints(o)
            wayPoints = o.K;
        end
        
        function K = getStrike(o)
            K = o.K;
        end
        
        function K = getContractSize(o)
            K = o.contractSize;
        end    
        
        function K = getCommission(o)
            K = o.commission;
        end          
        % Print out the instrument returning a string
        function name = print( instrument)        
            name = sprintf('Call K=%d',instrument.K);
        end
    end
    
end

