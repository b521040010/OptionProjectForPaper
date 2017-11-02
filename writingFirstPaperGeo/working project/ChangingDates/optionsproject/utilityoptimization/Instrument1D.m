classdef Instrument1D < Instrument
    %INSTRUMENT1D An instrument in 1-D problem so it's price only
    %   depends upon a single risk factor
    
    properties
    end
    
    methods
        function o = Instrument1D(bid,ask,bidSize,askSize)
            o = o@Instrument(bid,ask,bidSize,askSize);
        end       
        
    end

    methods (Abstract)
        delta = deltaAtInfinity(o)
        % As the stock tends to infinity, what is the delta? This is
        % used to ensure that liability is bounded
        wayPoints = getWaypoints(o)
        % Waypoints that should be included in the integration when
        % evaluating this option
    end
    
    
end

