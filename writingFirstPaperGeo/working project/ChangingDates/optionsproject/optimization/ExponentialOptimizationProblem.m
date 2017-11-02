classdef ExponentialOptimizationProblem < OptimizationProblem
    
    properties
        boundedLiability
        gamma
        
        % instruments that should be considered in the problem - those
        % which have no measurable value are excluded
        instrumentFilter
        
        % parts of the mosek optimization problem specification coming
        % from the objective function
        aObj
        cObj
        
        % parts of the mosek optimization problem specification coming
        % from the positivity constraints
        aConsPositive
        cConsPositive
        mapConsPositive
        
        % parts of the mosek optimization problem specification coming
        % from the cost constraint
        aConsCost
        cConsCost
        
        % An additional constraint given by requiring that we do not
        % take an unlimited liability
        aConsBoundedLiability
        cConsBoundedLiability
        
        % combined problem formulation
        geometricProgram
    end
    
    methods
        
        function o = ExponentialOptimizationProblem( dayData, model, boundedLiability )
            if (nargin<3)
                boundedLiability = false;
            end
            o = o@OptimizationProblem(dayData,model);
            o.boundedLiability = boundedLiability;
        end       
        
        function [ newPortfolio, utility ] = optimizeExponential(o, gamma)
            o = o.computeConstants(gamma);
            [newPortfolio,utility] = o.performOptimization();
        end
        
        function o = computeConstants( o, gamma )
            % We must write the problem in the form
            %  Minimize \sum_{k \in J_0} c^0_k \exp( a^0_k . x )
            % subject to
            %  \sum_{k \in J_i} c^i_k \exp( a^i_k . x ) < 1
            %
            % We let the vector x contain the quantities of each
            % instrument scaled by the absolute cost        
            %
            % By approximating the expected utility with a quadrature
            % rule we obtain our objective. So we get a values
            % c^0_k as k runs through the quadrature points S_k with
            % weights w_k. We wish
            % to maximize
            %    \int p(S) u( payoff if stock = S ) d S
            % =  \int p(S) u( x . payoff vector at S + payoff_0 ) d S
            % approx = \sum_k w_k p(s_k) u( x. payoff vec at S_k + payoff_0)
            %
            % We write b_k for the payoff vec when the stock price = S_k.
            % That is the ith component of b_k is the payoff of 1 dollars
            % worth of instrument i if the stock price = S_k compared with
            % the payoff of the initial portfolio, payoff_0.
            %
            % So we wish to maximize
            %   sum_k w_k p(s_k) u( x. b_k)
            % = sum_k w_k p(s_k) /gamma * (1- exp( -gamma (x. b_k + payoff_0) ))
            % Equivalently we wish to minimize
            %   sum_k w_k p(s_k)/gamma * exp( -(gamma x. b_k + payoff_0))
            % So we take a^0_k = -gamma b_k
            % and c^0_k = w_k p(s_k) exp( -gamma payoff_0)/gamma.
            %            
            
            o.gamma = gamma;
            o.utilityFunction = ExponentialUtilityFunction(gamma);
            quadRule = o.getQuadRule();
            
            % Use Mosek's special features to optimize an exponential utility
            % function
            probs =  o.model.pdf(quadRule.x);     
            probs = fixupRowVector('pdf',probs,quadRule.x);
            nPoints = length(quadRule.x );

            costs = o.ddd.instrumentPrice;
            nFiltered = length(costs);
            payoff0 = o.computePayoff(zeros(1,nFiltered),quadRule.x');

            o.aObj = zeros( nPoints, nFiltered );
            o.cObj = 1/gamma * (quadRule.weights' .* probs)' .* exp(-gamma * payoff0);

            for j=1:nFiltered
                q = zeros(1,nFiltered);
                q(j)=1;
                payoffs =  o.computePayoff(q,quadRule.x');
                payoffs = fixupRowVector('payoffs',payoffs - payoff0,quadRule.x) ;                
                o.aObj(:,j) = -gamma / abs(costs(j)) * payoffs';
%                fprintf( '%d Instrument strike=%d, type=%d, price=%d. Range1 (%d,%d) Range2 (%d,%d) \n', j, o.dayData.instrumentStrike(j), o.dayData.instrumentType(j), o.dayData.instrumentPrice(j), min(payoffs), max(payoffs), min(o.aObj(pointFilter,j)), max(o.aObj(pointFilter,j)) );
            end    
            
            % Some stock prices are so unlikely that cObj takes the value
            % 0. We wish to eliminate these points from the quad rule.
            % These filtered out points are stored in point filter            
            pointFilter = ~(o.cObj < 1e-8 | isnan(o.cObj));
            % Because some instruments only have value at filtered out
            % points, we need to filter out these instruments too            
            o.instrumentFilter = max(o.aObj(pointFilter,:),[],1)~=0 | min(o.aObj(pointFilter,:),[],1)~=0;
                        
            o.aObj=o.aObj(pointFilter,o.instrumentFilter);
            o.cObj=o.cObj(pointFilter);
                        
            nFiltered = sum( o.instrumentFilter );

            % The condition that the components of x are all positive can be written
            % \exp( -\delta^i . x ) < 1
            % So we get nInstruments constraints with
            % a^i = -\delta^i and c^i = 1.
            o.aConsPositive = -eye(nFiltered);
            o.cConsPositive = ones(nFiltered,1);
            o.mapConsPositive = (1:nFiltered)';
            
            % The condition that the total cost is less zero can be written as
            % the appropriately signed sum of the elements of x is less
            % than 0.
            % Let a_max = the sign of each instrument. Our condition is
            %     a_max .x < 0
            % Equivalently
            %     1 exp( a_max .x ) < 1
            % So take c_max = 1 and a_max as above.
            o.aConsCost = costs( o.instrumentFilter ) ./ abs(costs(o.instrumentFilter));
            o.cConsCost = 1;
            
            % The condition that we do not take on an unlimited liability
            unboundedInstruments = (o.ddd.instrumentType==DayData.callType | o.ddd.instrumentType==DayData.futureType);
            instrumentSigns = 2*(o.ddd.instrumentPrice>=0) -1;
            terminalSlopes = unboundedInstruments .* (instrumentSigns ./ abs(o.ddd.instrumentPrice));
            o.aConsBoundedLiability = -terminalSlopes( o.instrumentFilter );
            o.cConsBoundedLiability = 1;

            gp = GeometricProgram();
            gp = gp.setObjective( o.aObj, o.cObj );
            gp = gp.addConstraints( o.aConsPositive, o.cConsPositive, o.mapConsPositive );
            gp = gp.addConstraint( o.aConsCost, o.cConsCost );
            if o.boundedLiability
                gp = gp.addConstraint( o.aConsBoundedLiability, o.cConsBoundedLiability );
                gp = gp.addConstraints(eye(nFiltered), exp(-1)*ones(nFiltered,1), (1:nFiltered)');
            end  
            o.geometricProgram = gp;
            % Note that there is a trivial feasible point: x = 0.
        end
        
        function [utility, qScaled] = computeUtilityUsingMatrices(o, q )
            % For testing compute the utility using cObj and aObj            
            qScaled = q .*  abs(o.ddd.instrumentPrice);
            qDash = qScaled( o.instrumentFilter );
            expectedDisutilityMainTerm = sum(o.cObj .* exp(o.aObj * qDash'));
            utility = -expectedDisutilityMainTerm + 1/o.gamma;
        end               
        
        function validateObjective( o, index )
            costs = o.ddd.instrumentPrice;
            nInstruments = length(costs);
            testQ = zeros(1, nInstruments );
            testQ(index)=1;

            quadRule = o.getQuadRule();
            probs =  o.model.pdf(quadRule.x);     
            probs = fixupRowVector('pdf',probs,quadRule.x);
            
            unscaledQuantities = testQ ./ abs(costs);
            
            payoffs =  o.computePayoff(unscaledQuantities,quadRule.x');
            
            testQFiltered = testQ( o.instrumentFilter );
            disutility = exp(-o.gamma.*payoffs);
            expectedDisutilityMT = 1/o.gamma .* sum(quadRule.weights' .* probs .* disutility');
            expectedDisutilityMainTerm = sum(o.cObj .* exp(o.aObj * testQFiltered'));
            assertApproxEqual(expectedDisutilityMT, expectedDisutilityMainTerm, 1e-6);
            expectedUtility = o.computeExpectedUtilityDQ( unscaledQuantities);
            numericIntegral = o.computeExpectedUtilityDQ( unscaledQuantities, quadRule);
            expectedUtilityMainTerm = numericIntegral - 1/o.gamma;
            assertApproxEqual( expectedUtility, numericIntegral, 1e-4);
            assertApproxEqual( -expectedUtilityMainTerm, expectedDisutilityMainTerm, 1e-4 );
        end
        
        function dispMatrices(o)
            % Print out the various matrices for debugging purposes
            disp('a');
            formatMathematica(o.a);
            disp('c');
            formatMathematica(o.c);
            disp('map');
            formatMathematica(o.map);
        end
        
        function [newPortfolio, utility] = performOptimization(o)
            
            
            [~,res]  = mosekopt('symbcon'); 
            sc       = res.symbcon; 
            param = [];
            param.MSK_DPAR_INTPNT_TOL_DFEAS = 1e-3;
            param.MSK_IPAR_BI_IGNORE_MAX_ITER   = sc.MSK_ON; 
            param.MSK_IPAR_LOG = 0;
            gp = o.geometricProgram;
            res = mskgpopt(gp.c,gp.a,gp.map, param);
            %res = mskgpopt(gp.c,gp.a,gp.map);     
            if (~strcmp(res.sol.itr.solsta,'OPTIMAL') && ~strcmp(res.sol.itr.solsta,'NEAR_OPTIMAL'))
                % Repeat the optimization with debug information
                param.MSK_IPAR_LOG = 10;
                mskscopt(o.opr,o.opri,o.oprj,o.oprf,o.oprg,o.c,o.A,o.blc,o.buc,o.blx, o.bux, param);                                     
                error('Unable to find solution to optimization problem: %s',res.sol.itr.solsta);
            end
            costs = o.ddd.instrumentPrice;
            qFiltered = res.sol.itr.xx' ./ abs(costs(o.instrumentFilter)); 
            q = zeros(1,o.ddd.nInstruments);
            q(o.instrumentFilter) = qFiltered;
            utility = o.computeExpectedUtilityDQ( q );
            newPortfolio = o.constructNewPortfolio(q).simplify();
            
        end
                               
    end
    
    
end

