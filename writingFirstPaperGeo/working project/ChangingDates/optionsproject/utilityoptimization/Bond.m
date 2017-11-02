classdef Bond < Instrument1D
    %Bond Represents holding a risk free bond
    
    properties
        T
        r
    end

    properties (Constant)
        contractSize=1;
      commission=0;
    end
    
    methods
        function o = Bond(T, r, bid,ask,bidSize,askSize)
            o = o@Instrument1D(bid,ask,bidSize,askSize);
            o.ask=o.ask*o.contractSize;
            o.bid=o.bid*o.contractSize;
            o.T = T;
            o.r = r;
        end
        
        % Compute the payoff of an instrument in the given scenarios.
        % each row of the scenarios matrix represnts a specific scenario
        % so this method should return a column vector
        function value = payoff( o, scenarios ) 
            if o.bid>0
                price = o.bid;
            else 
                price = o.ask;
            end
            value = exp( o.r *o.T) * price * ones( size(scenarios));
        end
        
        function newBond = updateMaturity( bond, daysPassed )
            t = daysPassed/365;
            factor = exp( t * bond.r );
            newBond = Bond( bond.T-t, bond.r, factor*bond.bid, factor*bond.ask, bond.bidSize, bond.askSize );
        end
        
        function delta = deltaAtInfinity(~)
            delta = 0;
        end
        
        function wayPoints = getWaypoints(~)
            wayPoints = [];
        end
        
        
        function name = print( ~ )        
            name = 'Bond';
        end
        
    end
    
end

