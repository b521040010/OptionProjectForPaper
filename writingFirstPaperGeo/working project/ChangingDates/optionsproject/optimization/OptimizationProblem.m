classdef OptimizationProblem
    %OPTIMIZATIONPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dayData; % Price data today
        ddd; % Expanded day data
        model; % The beliefs about the future
        utilityFunction; % Utility to optimize
        currentPortfolio; % The current portfolio
        simpleInterest; % The interest rate
    end
    
    methods
        
        function o = OptimizationProblem( dayData, model )
            initutilityoptimization();
            o.dayData = dayData;
            o.ddd = DoubledDayData(dayData);
            o.model = model;
                        
            o.utilityFunction = ExponentialUtilityFunction(1);
            o.currentPortfolio = OldPortfolio();
            o.currentPortfolio.strikes = [];
            o.currentPortfolio.instrumentTypes = [];
            o.currentPortfolio.quantities = [];
            o = init(o);
        end       

        function o = init( o)
            interest = o.dayData.getInterestRate();
            o.simpleInterest = exp(o.model.T * interest) - 1.0;  
        end
        
        
        function newPortfolio = constructNewPortfolio(o,q)
            % The portfolio that will result from adding the given
            % quantities
            newPortfolio = OldPortfolio();
            newPortfolio.strikes = [ o.currentPortfolio.strikes o.ddd.instrumentStrike ];
            newPortfolio.instrumentTypes = [ o.currentPortfolio.instrumentTypes o.ddd.instrumentType ];
            sgn = 2.*(-o.ddd.isShort>=0) - 1;
            %sgn = 2.*(o.ddd.instrumentPrice>0) - 1;
            newPortfolio.quantities = [ o.currentPortfolio.quantities (q .* sgn)];
        end        
        
        function r = computePayoff(o,q, scenarioAsks)
        % Given the quantities bought and sold computes
        % the resultant payoff
            newPortfolio = o.constructNewPortfolio(q);
            r = newPortfolio.payoff( o.simpleInterest, scenarioAsks );
        end        
                
        
        function [r] = computeExpectedUtility(o, netQ, quadRule )
            dq = o.ddd.doubleQ(netQ');
            if (nargin==3)
                r = computeExpectedUtilityDQ(o, dq, quadRule );
            else
                r = computeExpectedUtilityDQ(o, dq);
            end
        end
        
        function [ r ] = computeExpectedUtilityDQ( o, quantities, quadRule)
            % Compute the expected utility given the doubled quantities
            % column vector
            function r = integrand(ask)
                payoffs =  o.computePayoff(quantities,ask);                 
                
                logProbs = o.model.logPdf(ask);
                logProbs = fixupRowVector('logPdf',logProbs,ask);
                payoffs = fixupRowVector('netValue',payoffs,ask);
                r = o.utilityFunction.weightedUtility(payoffs, logProbs);  
                if (sum(isfinite(r))~=length(ask)) 
                    disp('quantities');
                    disp(quantities);
                    disp('payoffs');
                    disp(payoffs);
                    disp('logProbs');
                    disp(logProbs);
                    disp('Weighted utility');
                    disp(r);
                    disp('Ask');
                    disp(ask);
                    error('NaN encountered');
                end
                    
            end

            % Numeric intergration does not perform well without some clues as
            % to the location of the probability mass
            
            if nargin<3
                wayPoints = [ o.model.getWayPoints() o.ddd.instrumentStrike o.currentPortfolio.strikes];
                wayPoints = sort( unique(wayPoints ));
                s = integral(@integrand, 0, Inf, 'WayPoints', wayPoints);
            else
                s = quadRule.integrate(@integrand);
            end
            r=s;

        end        
        
        function [ newPortfolio, utility, q ] = optimizeFmincon(o, netQ0)
        %OPTIMIZEUTILITYGENERIC Performs utility optimization given the
        % payoff function, utility function and pdf. Optimization performed
        % using fmincon. You can suggest some starting quantities if you
        % like
        
            quadRule = o.getQuadRule();
        
            function r = disutilityForQuantities(q)
                r = -o.computeExpectedUtilityDQ( q, quadRule );
            end

            costVector = o.ddd.instrumentPrice;
            nInstruments = length( costVector );
            q0 = zeros(1,nInstruments);
            if (nargin>=2)
                q0 = o.ddd.doubleQ( netQ0' );
            end
            A = [ -eye(nInstruments) ; costVector ];
            b = [ zeros(1,nInstruments) 0 ];            
            options = optimset('fmincon');
            options = optimset(options,'Display','off');
            options = optimset(options,'Algorithm', 'active-set');
            
            %options = optimset(options,'TolFun',1e-3);
            [q, ~, exitFlag] = fmincon( @disutilityForQuantities, q0, A, b, [],[],[],[],[], options );
            if ( exitFlag < 0 )
                options = optimset(options,'Display','on');
                fmincon( @disutilityForQuantities, q0, A, b, [],[],[],[],[], options );
                error('fmincon failed');
            end
            q = q .* (q > 0); % Get rid of small errors that result in some negative q values

            utility = -disutilityForQuantities(q);
            newPortfolio = o.constructNewPortfolio(q).simplify();
            
            q = o.ddd.netQ(q)';
        end        
                        
        
        function quadRule = getQuadRule(o)
            % Get an accurate quadrature rule adapted to this problem
            quadRule = QuadRule.adapted( @(x) o.model.pdf(x), o.ddd.instrumentStrike, o.currentPortfolio.strikes, o.model.getWayPoints() );
        end
        
        function ump = createUtilityMaximizationProblem(op)
            % Create a utility maximization problem corresponding
            % to this problem
            ump = UtilityMaximizationProblem1D();
            ump.setModel( op.model );
            
            ump.setUtilityFunction(op.utilityFunction);

            dd = op.dayData;
                        
            for i=1:length(dd.instruments)
                ump.addInstrument( dd.instruments{i} );
            end
            
            pAsI = PortfolioAsInstrument( op.currentPortfolio, op.dayData, op.simpleInterest );
            ump.setCurrentPosition( pAsI );
        end
        
                               
    end
    
    
end

