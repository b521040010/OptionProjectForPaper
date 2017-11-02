% A rescaled version of another distribution with given
% mean and standard deviation
classdef rescaledDistribution
    properties
        mu
        sigma
        delegate
    end
    methods
        % Integrate (at + b)*pdf from t1 to t2
        function res = integrate( obj, a, b, t1, t2)
            s1 = (t1 - obj.mu)/obj.sigma;
            s2 = (t2 - obj.mu)/obj.sigma;
            res = obj.delegate.integrate( a*obj.sigma, a*obj.mu + b, s1,s2 );
        end
    end
end
