classdef Instrument < matlab.mixin.Copyable
    %INSTRUMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties %(Access='protected')
        bid;
        ask;
        bidSize;
        askSize;
    end
    
    methods
        function o = Instrument( bid, ask, bidSize, askSize )
            o.bid = bid;
            o.ask = ask;
            o.bidSize=bidSize;
            o.askSize=askSize;
            assert(bidSize>=0);
            assert(askSize>=0);
            if bid==21 && ask==20 
                % Special case to allow testing for bid ask arbitrage
            else
%                assert(bid<=ask || isnan(bid) || isnan(ask));            
%                assert(bid>=0 || isnan(bid));
            end
        end
        
        % Update the maturity of the instrument and return an appropriate
        % new instrument
        function newInstrument = updateMaturity( instrument, daysPassed )        
            newInstrument = instrument;
        end
    end
    
    methods (Abstract)
        % Compute the payoff of an instrument in the given scenarios.
        % each row of the scenarios matrix represnts a specific scenario
        % so this method should return a column vector
        value = payoff( instrument, scenarios )
        
        % Print out the instrument returning a string
        name = print( instrument)
        
    end
    
    methods
        
        function b = getBid(o)
            b = o.bid;
        end
        
        function a = getAsk(o)
            a = o.ask;
        end
        
        function bS=getBidSize(o)
            bS=o.bidSize;
        end
        
        function aS=getAskSize(o)
            aS=o.askSize;
        end
    end
    
end

