classdef PowerUtilityFunction < UtilityFunction
    %POWERUTILITYFUNCTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        eta % Power utility parameter
        name % Component understood by Mosek
        f    % Component understood by Mosek
        g    % Component understood by Mosek
        c    % Component understood by Mosek
    end
    
    methods
        function o = PowerUtilityFunction(eta)
            o.eta = eta;
            if (o.eta==1)
                o.name = 'log';
                o.f = 1;
                o.g = 1;
                o.c = 0;
            else
                o.name = 'pow';
                o.g = 1 - o.eta;
                o.f = 1/(1-o.eta);
                o.c = -o.f;
            end
        end
        
        % Gets the name to be passed to Mosek for separable optimization
        function name = getMosekName(o)
            name = o.name;
        end
        % Get the value of f to pass to Mosek for separable optimization
        function f = getF(o, logProbs)
            f = o.f .* exp(logProbs);
        end
        % Get the value of g to pass to Mosek for separable optimization
        function g = getG(o)
            g = o.g;
        end 
        % Get the constant to add onto the utility as understood by Mosek
        function c = getConstant(o)
            c = o.c;
        end
        % Get the lower bound on the payoff for this utility function
        function lb = getLowerBound(o)
            lb = 0;
        end
        
        function r = weightedUtility( o, x, logProb) 
            if (o.eta==1)
                r = log(x) .* exp( logProb );
            else
                r = (x.^(1-o.eta) - 1)/(1-o.eta) .* exp(logProb);
            end                
            r(x<0)=-Inf;
        end
            
        
    end
    
end

