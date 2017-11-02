classdef UtilityFunction
    %UTILITYFUNCTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Abstract)
        % Gets the name to be passed to Mosek for separable optimization
        name = getMosekName(o)
        % Get the value of f to pass to Mosek for separable optimization
        % given the log probability.  This method
        % should be used in conjuction with adjustXj.
        % The mosek value computed from xj = x+adjustXj( logProb)
        %                               f = getF(o, logProb)
        %                               g = getF(o)
        % will be called mVal.
        % Should equal be e^logProb (u(x)) = mVal + e^logProb * getConstant();
        % where u is the utility function.
        f = getF(o, logProb)
        % Get the value of g to pass to Mosek for separable optimization
        g = getG(o)
        % Get the constant to add onto the utility as understood by Mosek
        c = getConstant(o)
        % Get the lower bound on the payoff for this utility function
        lb = getLowerBound(o);
        % Compute the utility at x times e^logProb. The strange signature
        % ensures better numerical behaviour for exponential utility
        r = weightedUtility( o, x, logProb);
    end
    
    methods
        function testEvaluation( o, x )
            name = o.getMosekName();
            logProb = log(0.7)*ones(size(x));
            adjustment = o.xjAdjustment( logProb );
            xj = x+adjustment;
            f = o.getF(logProb);
            g = o.getG();
            c = getConstant(o);
            val1 = SeparableProblem.sepEval( name, xj,f ,g ) + exp(logProb)*c;
            val2 = weightedUtility( o, x, logProb);
            assertApproxEqual( val1, val2, 0.0001);
        end
        
        function adjustment = xjAdjustment( o, logProb )
            % For better numerical behaviour it is a good idea to adjust
            % the payoff variables to automatically include the probability
            % multiplier
            adjustment = zeros(size(logProb));
        end
    end
    
end

