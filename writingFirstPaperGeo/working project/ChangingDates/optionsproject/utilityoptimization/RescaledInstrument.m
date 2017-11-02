classdef RescaledInstrument < Instrument
    
    properties 
        scale
        delegate
    end        
    
    methods
        
        function o = RescaledInstrument( instrument, scale )
            o@Instrument(scale*instrument.getBid(),scale*instrument.getAsk(),instrument.getBidSize(),instrument.getAskSize());
            o.delegate = instrument;
            o.scale = scale;
        end
        
        function delta = deltaAtInfinity(o)
            delta = o.scale * o.delegate.deltaAtInfinity();
        end   
        
        % Compute the payoff of an instrument in the given scenarios.
        % each row of the scenarios matrix represnts a specific scenario
        % so this method should return a column vector
        function value = payoff( o, scenarios )
            value = o.scale * o.delegate.payoff( scenarios );
        end
        
        function name = print( o )        
            name = sprintf('%d * %s', o.scale,o.delegate.print());
        end        
          
    end
end
