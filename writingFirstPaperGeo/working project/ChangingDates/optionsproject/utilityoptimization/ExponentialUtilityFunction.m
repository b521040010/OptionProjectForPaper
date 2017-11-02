classdef ExponentialUtilityFunction < UtilityFunction

    properties
        gamma
    end

    methods
        
        function o = ExponentialUtilityFunction(gamma)
            o.gamma = gamma;
        end        
        
        function [ r ] = weightedUtility( o, x, logProb) 
            % Compute the utility weighted with the given log probabilities
            if (nargin==2)
                logProb = 0;
            end
            r = -1/o.gamma*(exp( - o.gamma * x + logProb) - exp( logProb ));
        end
                
        % Gets the name to be passed to Mosek for separable optimization
        function name = getMosekName(o)
            name = 'exp';
        end
        % Get the value of f to pass to Mosek for separable optimization
        function f = getF(o, logProbs)
            %f = -1;
            f = -1/o.gamma;
        end
        % Get the value of g to pass to Mosek for separable optimization
        function g = getG(o)
            g = -o.gamma;
        end 
        % Get the constant to add onto the utility as understood by Mosek
        function c = getConstant(o)
            %c = 1;
            c = 1/o.gamma;
        end
        % Get the lower bound on the payoff for this utility function
        function lb = getLowerBound(o)
            lb = -Inf;
        end
        
        function adjustment = xjAdjustment( o, logProb )
            % For better numerical behaviour it is a good idea to adjust
            % the payoff variables to automatically include the probability
            % multiplier
            adjustment = -1/o.gamma * logProb;
        end        
        
    end

end

