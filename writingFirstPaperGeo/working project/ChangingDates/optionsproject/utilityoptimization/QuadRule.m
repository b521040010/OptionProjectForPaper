classdef QuadRule
    % A quadrule is a quadrature rule over a certain interval
    properties
        a;
        b;
        weights;
        x;
    end

    methods( Static)
        function quadRule = combine( quadRules )
            % Compose a sequence of consecutive quadrature rules
            % into one large rule
            quadRule = QuadRule();
            quadRule.a = quadRules(1).a;
            quadRule.b = quadRules(1).b;
            totalLength = 0;
            for i=1:length(quadRules)
                qr = quadRules(i);
                quadRule.a = min( quadRule.a, qr.a );
                quadRule.b = max( quadRule.b, qr.b );
                totalLength = totalLength + length( qr.weights );
            end
            quadRule.weights = zeros(totalLength,1);
            quadRule.x = zeros(totalLength,1);
            current = 1;
            for i=1:length(quadRules)
                qr = quadRules(i);
                len = length( qr.weights );
                quadRule.weights( current:current+len-1) = qr.weights;
                quadRule.x( current:current+len-1) = qr.x;
                current = current + len;
            end
        end
        
        
        function quadRule = gaussLegendre(n, a, b)
            if nargin==1
                quadRule = QuadRule();
                quadRule.a = -1;
                quadRule.b = 1;

                % Construct the legendre polynomial
                P = legendrePolynomial(n);
                quadRule.x = sort(roots( P(end:-1:1)));
                assert( length(quadRule.x)==n);
                quadRule.weights = zeros(n,1);

                p1 = sum( P );

                PDash = P(2:end) .* (1:n);
                for k=1:n
                    xk = quadRule.x(k);
                    xVec = xk .^ (0:n-1);
                    pd = sum(xVec .* PDash);
                    quadRule.weights(k) = 2/((1-xk^2)*((pd/p1)^2));
                end
            else
                quadRule = QuadRule.gaussLegendre(n);
                quadRule = quadRule.transform(a,b);
            end
        end       
                
        function quadRule = adapted(varargin) 
            % Create a quad rule adapter to the given points
            points = [];
            pdf = varargin{1};
            for i=2:nargin
                points = [points varargin{i}];
            end
            points = sort( unique( points ));
            start = points(1);
            for i=1:(length(points)-1)
                quadRules(i) = QuadRule.gaussLegendre(5,start,points(i+1));
                %val = quadRules(i).integrate(pdf);
                %if (val<1e-6)
                %    quadRules(i)=QuadRule.gaussLegendre(1,start,points(i+1));
                %end
                start = points(i+1);
            end
            quadRule = QuadRule.combine( quadRules );
        end        

    end        
    
    methods
        
        function ret = integrate( quadRule, f )
            % Integrate a function using the quadrature rule
            ret = sum(f(quadRule.x') .* quadRule.weights');
        end
        
        function quadRule = transform(quadRuleIn, a, b) 
            quadRule = QuadRule();
            quadRule.a = a;
            quadRule.b = b;
            quadRule.x = (quadRuleIn.x-quadRuleIn.a)*(b-a)/(quadRuleIn.b-quadRuleIn.a) + a;
            quadRule.weights = quadRuleIn.weights*(b-a)/(quadRuleIn.b-quadRuleIn.a);
        end
        
    end
    
end

