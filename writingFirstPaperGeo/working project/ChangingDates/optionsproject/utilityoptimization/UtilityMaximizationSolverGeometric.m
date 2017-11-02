classdef UtilityMaximizationSolverGeometric < handle
    %UTILITYMAXIMIZATIONSOLVER This class builds up a separable problem
    %  matching a utility maximization problem.
    %  The variables are
    %    q  a vector of positive quantities of each instrument
    %    P  the payoffs at each x point in the quadrule of the problem
    %    xj  the value sent to Mosek's separable optimization computed
    %        as P with a possible offset to improve numerics for
    %        exponential optimization
    %        By storing the variables in this neat fashion, it is easy
    %        to add constraints on the payoffs are the quantities. 
    %
    %  There are quite a few issues one needs to consider to get this
    %  to work.
    %  
    %  Firstly equality constraints don't appear to work very
    %  well, so wherever possible (everywhere) a single equality is chosen.
    %
    %  Secondly, this is nonlinear optimization so scaling is important.
    %  The rescaling is performed in UtilityMaximizationProblem before the
    %  solver is called to simplify the solver class.
    %
    %  Thirdly, when comparing with fmincon one has to accept that fmincon
    %  is imperfect and can give terrible answers. How depressing!
    %
    %  Fourthly, when performing exponential utility optimization it is
    %  good to incorporate the probabilities into the xj variables
    %
    %  Fifthly, adding a condition that the maximum liability is bounded
    %  improves exponential utility optimization
    %
    %  Sixth, one must add reasonable tolerance parameters to the solver.
    %
    %  Seventh, instruments that don't have a reasonable probability of a
    %  profit/loss result in an ill posed problem, so we must exclude them
    %  from consideration
    %
    %  The code is made much easier to follow by some good structuring. We
    %  allow constraints to be added in separate pieces rather than a
    %  single giant constraint. We also allow them to be labelled so one
    %  can test if constraints pass and get helpful messages about the
    %  source of the problem.
    
    
    properties
        % The problem we are solving        
        p 
        % The scenarios we will actually consider - note this is a hangover
        % from a failed version of the code - it is important to consider
        % all scenarios
        scenarioFilter
        % The instrument corresponding to the ith investment choice and
        % whether this is a short position. Only instruments that have
        % non zero payoffs will be considered
        indexToInstrument
        indexToShort
        % A cutoff to use to determine which filters need to be considered
        logProbCutoff
        % The number of quantity variables
        nq
        % The number of payoff variables
        nP
        % Mosek separable program object
        geoProblem
        % Payoff matrix. When q is multiplied by this we get the payoffs
        % modulo offset p.payoff0
        payoff
        % The vector of costs
        costVec
    end
    
    methods
        
        function o = UtilityMaximizationSolverGeometric( p )
            o.p = p;
            o.logProbCutoff = log(1e-9);
     %       o.logProbCutoff = -1;
            o.selectScenariosAndInstruments();
            %o.computeObjectiveFunction();
            o.addCostConstraint();
            o.addPqIdentityConstraint();
            %o.addXjConstraint();
            o.addPositiveQuantityConstraint();
            nConstraints = length(p.constraints);
            for i=1:nConstraints
                constraint = p.constraints{i};
                constraint.applyConstraint(o, o.geoProblem );
            end
        end
        
       
        
        function o = selectScenariosAndInstruments(o)
            o.scenarioFilter = o.p.logProb>o.logProbCutoff;
            
            nInstruments = length(o.p.instruments);
            count = 1;
            for i=1:nInstruments                
                payoffi = o.p.payoffi{i};
                % Ignore worthless instruments
                if max(abs(payoffi(o.scenarioFilter)))>0
                    o.indexToInstrument(count) = i;
                    o.indexToShort(count) = true;
                    o.indexToInstrument(count+1) = i;
                    o.indexToShort(count+1) = false;
                    count = count+2;
                end
            end
            
            % filtering out scenarios causes problems...
            o.scenarioFilter = o.p.logProb>-Inf;
            o.nq = length(o.indexToInstrument);
            o.nP = sum( o.scenarioFilter);
            o.geoProblem = GeometricProblem(o.nq+1);
            
        end

        function addPqIdentityConstraint(o)
            % The variables q and P are related by the payoffs.
            % This adds a constraint ensuring that relation holds (or
            % rather that the payoff is no bigger than the max
            o.payoff = zeros(o.nP,o.nq);
            for i=1:o.nq
                instrument = o.indexToInstrument(i);
                short = o.indexToShort(i);
                sgn = 1-2*short;
                o.payoff(:,i)=o.p.payoffi{instrument}(o.scenarioFilter) * sgn;
                if ~short
                    totalPayoff = o.payoff(:,i)+o.payoff(:,i-1);
                    assert( sum( totalPayoff>0 )==0 );
                end
            end       
            %payy=sparse(horzcat(o.payoff,o.p.payoff0(o.scenarioFilter)));
            payy=horzcat(o.payoff,o.p.payoff0(o.scenarioFilter));
            logProb = o.p.logProb(o.scenarioFilter);
            size(payy);
            logProb;
            payy=sparse(payy);
            o.geoProblem.setObjective(exp(logProb),(-o.p.utilityFunction.gamma*payy));
            o.geoProblem.addLinearConstraint( horzcat(zeros(1,o.nq),1) , 1, 'to make the last variable =1 , for the payoff0' ); 
%             o.geoProblem.addLinearUpperBound(horzcat(zeros(1,o.nq),1) , 1+0.00001,'cost constraint');
%             o.geoProblem.addLinearLowerBound(horzcat(zeros(1,o.nq),1) , 1-0.00001,'cost constraint');
        end
        
        function addCostConstraint(o)
            % The total cost of the portfolio must be OK
            vec = zeros(1,o.nq+1);
            for i=1:o.nq
                instrumentIdx = o.indexToInstrument(i);
                instrument = o.p.instruments{instrumentIdx};
                if (o.indexToShort(i))
                    vec(i) = -instrument.getBid();
                else
                    vec(i) = instrument.getAsk();
                end
                % If the price is infinite we must add the constraint
                % that the instrument cannot be bought/sold at all.
%                 if (~isfinite(vec(i))) 
%                     vec(i) = 0;
%                     o.sepProblem.bux(i)=0;
%                 end
            end
            o.costVec = vec(1:o.nq);
            o.geoProblem.addLinearUpperBound(horzcat(o.costVec,0),0,'cost constraint');

        end    
        
        function addPositiveQuantityConstraint(o)
            
            A = horzcat(eye(o.nq), zeros(o.nq,1));
            b = zeros(o.nq,1);
            o.geoProblem.addLinearLowerBound(A,b, 'Positive quantity');
        end            
        
%         function computeObjectiveFunction(o)
%             
%             logProb = o.p.logProb(o.scenarioFilter);
%             % Computes the objective function
%             name = o.p.utilityFunction.getMosekName();
%             f = o.p.utilityFunction.getF( logProb );
%             g = o.p.utilityFunction.getG();
%             
%             c = zeros( o.nq + o.nP, 1 );
%             opr = repmat( name, o.nP, 1 );
%             oprj = o.nq  + (1:o.nP)';
%             oprf = -f.*ones(o.nP,1); % We minimize disutility, hence the minus sign
%             oprg = g*ones( o.nP, 1 );
%             o.sepProblem.setObjective(c, opr,oprj,oprf,oprg);
%         end

%         function utility = utilityForQuantities(o, netQ)            
%             % Compute the utility for the given quantities
%             x = o.computeX(netQ);
%             objValue = o.sepProblem.computeObjective(x');
%             utility = -objValue + o.p.utilityFunction.getConstant();
%         end
%         
%         function assertConstraintsPassed(o, netQ, tolerance) 
%             % Assert all constraints pass
%             x = o.computeX(netQ);
%             o.sepProblem.assertConstraintsPassed(x, tolerance);
%         end
%         
%         
%         function x = computeX(o, netQ )
%             % Compute the vector of q + payoffs given
%             % the net value of q
%             assert( size(netQ,2)==1); % should be a column vector
%             q = zeros( o.nq, 1);
%             for idx=1:o.nq
%                 instrument = o.indexToInstrument(idx);
%                 isShort = o.indexToShort(idx);
%                 if (isShort && netQ(instrument)<0)
%                     q(idx)=-netQ(instrument);
%                 elseif (~isShort && netQ(instrument)>0)
%                     q(idx)=netQ(instrument);
%                 end
%             end
%             temp=o.p.currentPosition.map;
%             temp=values(temp);
%             r=temp{1}.instrument.r;
%             T=temp{1}.instrument.T;
%             S0erT=o.costVec*exp(r*T);
%             S0erT=repmat(S0erT,o.nP,1);
%             P = o.payoff * q + o.p.payoff0(o.scenarioFilter)-S0erT*q;
%             %P = o.payoff * q + o.p.payoff0(o.scenarioFilter);
%             
%             logProb = o.p.logProb(o.scenarioFilter);
%             adjustment = o.p.utilityFunction.xjAdjustment(logProb);
%             xj = P + adjustment;
%             x = vertcat(q,P, xj);  
%         end
        
%         function writeProblem( o, file )
%             o.sepProblem.writeProblem( file );
%         end

        function [utility,quantities,qp] = solve(o)
            % Solve the problem computing the resulting utility and 
            % the quantities to hold of each instrument                        
            o.geoProblem.tolerance = o.p.tolerance;
            o.geoProblem;
            size(o.geoProblem.w);
            min(o.geoProblem.w);
            size(o.geoProblem.u);
            o.geoProblem.u;
            size(o.geoProblem.ao);
            o.geoProblem.ao;
            size(o.geoProblem.ac);
            o.geoProblem.ac;
            size(o.geoProblem.c);
            o.geoProblem.c;
            [objective,qp] = o.geoProblem.optimize();  
            
            o.geoProblem.w;
            o.geoProblem.ao;
            o.geoProblem.ao*qp;
            exp(o.geoProblem.ao*qp);
            sum(o.geoProblem.w.*exp(o.geoProblem.ao*qp),1);
            log(sum(o.geoProblem.w.*exp(o.geoProblem.ao*qp),1));
            o.geoProblem.computeObjective(qp);
            o.geoProblem.assertConstraintsPassed(qp,0.001);
            %utility=objective;
            utility =(1-exp(objective))*o.p.utilityFunction.getConstant();
            nInstruments = size(o.p.instruments,2);
            quantities = zeros(nInstruments,1);
            for i=1:o.nq
                ins = o.indexToInstrument(i);
                short = o.indexToShort(i);
                sign = 1-2*short;
                quantities(ins) = quantities(ins)+ sign*qp(i);
            end
        end
        
    end
    
end


