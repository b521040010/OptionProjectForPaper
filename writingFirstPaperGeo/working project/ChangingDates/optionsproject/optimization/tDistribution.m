% Implementation of the student t distribution
classdef tDistribution
    properties
        nu
    end
    methods
        function obj = tDistribution()
        end
        % Integrate a polynomial from t1 to t2
        function res = integrate( obj, poly, t1, t2 )
            res = 0;
            for i=1:length(poly)
                res = obj.integratePower(i-1,t2)-obj.integratePower(i-1,t1);
            end
        end         
        function res = integratePower( obj, lambda, t ) 
            eta = obj.nu;
            res = (t^(1 + lambda)*(eta/(eta + t^2))^((1 + eta)/2)*((eta + t^2)/eta)^(( 1 + eta)/2) * hypergeom([(1 + eta)/2, (1 + lambda)/2], (3 + lambda)/2, -(t^2/eta)))/(sqrt(eta)*(1 + lambda)*beta(eta/2, 1/2));
        end            
    end
end
