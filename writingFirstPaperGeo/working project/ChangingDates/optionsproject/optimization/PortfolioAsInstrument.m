classdef PortfolioAsInstrument < Instrument
    %PORTFOLIOASINSTRUMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        portfolio
        simpleInterest
    end
    
    methods
        
        function pAsI = PortfolioAsInstrument( portfolio, dayData, simpleInterest )
            pAsI@Instrument(portfolio.value(dayData),portfolio.cost(dayData),0,0);
            assert( isfinite( pAsI.getAsk() ));
            assert( isfinite( pAsI.getBid() ));
            pAsI.portfolio = portfolio;
            pAsI.simpleInterest = simpleInterest;
        end
        
        % Compute the payoff of an instrument in the given scenarios.
        % each row of the scenarios matrix represnts a specific scenario
        % so this method should return a column vector
        function value = payoff( i, scenarios )
            value = i.portfolio.payoff( i.simpleInterest, scenarios' );
        end
        
        % Print out the instrument returning a string
        function name = print( instrument)
            name = 'Portfolio';
        end
        
        function waypoints = getWaypoints(o)
            waypoints = o.portfolio.strikes;
        end
    end
    
end

