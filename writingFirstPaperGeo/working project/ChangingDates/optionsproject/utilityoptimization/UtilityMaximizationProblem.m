classdef UtilityMaximizationProblem < matlab.mixin.Copyable
    % UtilityMaximizationProblem solves a utility maximization problem
    %   defined in terms of instruments
    %
    
    properties           
        
        instruments % The instruments available, stored as cell array
        % The current position represented as an instrument object
        currentPosition
        % The scenarios and weights are determined by
        % the quadrature rule and belief about future prices. Each
        % row of the scenario matrix represnts a scenario. The weight
        % is the associated "probability"
        scenarios
        logProb
        
        % The payoff of the current position
        payoff0
        % The payoff of the ith investible. Stored as cell array
        payoffi        
        
        % The utility function
        utilityFunction
        
        % Additional constraints
        constraints;
        
        % The duality gap tolerance for utility computations
        tolerance;
    end
    
    methods
        
        function o = UtilityMaximizationProblem()
            o.utilityFunction = PowerUtilityFunction(1);
            o.currentPosition = Bond(0,0,0,0,Inf,Inf);
            o.tolerance = 1e-9;
        end
        
        function addConstraint( o, c )
            idx = length( o.constraints )+1;
            o.constraints{idx} = c;
        end
        
        function setQuadRule( p, scenarios, logProb )
            % Sets the quadrature ruled used to solve the problem using
            % column vectors. For each scenario we have an associated log
            % probability of it occuring
            assert( size( scenarios,2 )==1);
            assert( size( logProb,2 )==1);
            assert( size( scenarios,1 )==size( logProb,1 ));
            p.scenarios = scenarios;
            p.logProb = logProb;
            p.payoff0 = p.currentPosition.payoff( p.scenarios );
            n = length( p.instruments );
            p.payoffi = cell(n,1);
            for i=1:n
                p.payoffi{i} = p.instruments{i}.payoff( scenarios );
            end
            
        end
        
        function p = indifferencePriceExponentialUtility( o, cashInstrument, quantity, instrument, priceGuess )
            % Perform buyers indifference pricing in the special case of
            % exponential utility and a totally liquid risk free bond
            assert( isa(cashInstrument,'Bond'));
            assert( isa(o.utilityFunction,'ExponentialUtilityFunction'));
            assert( cashInstrument.getBid()==cashInstrument.getAsk());
            
            V0 = o.optimize();
            copyProblem = copy(o);
            copyProblem.addToCurrentPosition( quantity, instrument);
            copyProblem.addToCurrentPosition( -priceGuess, cashInstrument );
            V1 = copyProblem.optimize();
            
            a = o.utilityFunction.gamma;
            p = exp(-cashInstrument.r*cashInstrument.T)*(1/a)*log((1-a*V0)/(1-a*V1)) + priceGuess;
        end

        
        function [p,quantities] = indifferencePrice( o, cashInstrument, quantity, instrument, guess )
            % Compute the buyers indifference price for k quantities of the instrument
            % This may be negative in order to compute bid prices
            V0 = o.optimize();
            copyProblem = copy(o);
            %we have only cash(initial wealth) for now
%             copyProblem.currentPosition.map
%             keys(copyProblem.currentPosition.map)
%            
%             temp=values(copyProblem.currentPosition.map)
%             temp{1}
            
            
            copyProblem.addToCurrentPosition( quantity, instrument);
%             disp('======')
%            copyProblem.currentPosition.map
            %now, we have calloption, portfolio of initial wealth
%             key=keys(copyProblem.currentPosition.map)
%             key{1}
%             key{2}
%             
%             value=values(copyProblem.currentPosition.map)
%             value{1}
%             value{2}
%             value{2}.instrument
%             value{2}.instrument.map
%             keys(value{2}.instrument.map)
%             values(value{2}.instrument.map)
%             temp2=values(value{2}.instrument.map)
%             temp2{1}
%             disp('==========')

            iterCount = 1;
            function [diff,quantities]=toSolve( priceGuess )
                resetData = copyProblem.addToCurrentPosition( -priceGuess, cashInstrument );
%                 disp('---------------------------------------------')
%                 copyProblem.currentPosition
%                 copyProblem.currentPosition.map
%                 key=keys(copyProblem.currentPosition.map)
%                 value=values(copyProblem.currentPosition.map)
%                 value{1}
%                 value{2}
%                 value{2}.instrument
%                 value{2}.instrument.map
%                 keys(value{2}.instrument.map)
%                 value=values(value{2}.instrument.map)
%                 value{1}
%                 temp2=values(value{2}.instrument.map)
%                 temp2{1}
%                 disp('---------------------------------------------')
                [V,quantities] = copyProblem.optimize();
                copyProblem.resetPosition( resetData );
                diff = V0-V;
                fprintf('Iteration %d\n', iterCount);
                iterCount = iterCount+1;
            end
            
            % Surprisingly bisection works better than using fzero
            % or fsolve (which give the same result)
            
            %options = optimset('fzero');
            %options = optimset(options,'Display','off');
            %options = optimset(options,'TolX', o.tolerance);
            %[negP,~,exitFlag] = fzero( @toSolve, -guess, options );
            %if (  exitFlag <= 0)
                % Smooth search failed try bisection
                %fprintf('Smooth search failed, using bisection\n');
                [p,quantities] = lineSearch( @toSolve, guess, o.tolerance);
            %else
            %    p = negP;
            %end
        end                
        
        function setCurrentPosition( p, currentPosition ) 
            % Sets the current position
            p.currentPosition = currentPosition;
            p.payoff0 = currentPosition.payoff( p.scenarios );
        end
        
        function resetData = addToCurrentPosition( p, quantity, instrument )
               % Add an instrument to the current position. Returns a
            % resetData object so this change can be reversed easily

            resetData.currentPosition = p.currentPosition;
            resetData.payoff0 = p.payoff0;

            port = Portfolio();
            additionalQ = ones(2,1);
            additionalQ(2) = quantity;
            additionalInstruments = cell(1,2);
            additionalInstruments{1}=p.currentPosition;
            additionalInstruments{2}=instrument;
            port.add(additionalQ, additionalInstruments);
            p.currentPosition = port;         
            p.payoff0 = p.payoff0 + quantity*instrument.payoff( p.scenarios ); 

        end
        
        function resetPosition( p, resetData )
            % Reverse a change made by addToCurrentPosition
            p.currentPosition = resetData.currentPosition;
            p.payoff0 = resetData.payoff0;
        end
        
        
        function addInstrument( p, i )
            % Adds a possible instrument
            idx = length( p.instruments )+1;
            p.instruments{ idx} = i;
            p.payoffi{idx} = i.payoff( p.scenarios );    
        end
        
        function removeInstrument( p, filter )
            error('Not implemented');
        end
        
        function setUtilityFunction( o, utilityFunction )
            o.utilityFunction = utilityFunction;
        end
        
        function setModel( o, model ) 
            % Convenience function to set the quad rule using a model
            % Get an accurate quadrature rule adapted to this problem
            nInstruments = length(o.instruments);
            wayPoints = [];
            for i=1:nInstruments
                instrument = o.instruments{i};
                wayPoints = horzcat(wayPoints,instrument.getWayPoints());
            end
            wayPoints = horzcat( wayPoints, o.currentPosition.getWayPoints());
            quadRule = QuadRule.adapted( wayPoints, o.model.getWayPoints() );
            x = quadRule.x;
            lp = log(weight) + model.logProb(wayPoints)

%             prices = model.simulatePricePaths(10000,1);
%             scenarios = prices(:,end);
%             wayPoints = scenarios';
%             x = wayPoints;
%             weight = (1/length(wayPoints))*ones(1,length(wayPoints));
%             lp = log(weight)
          
            
            o.setQuadRule( x, lp );
        end

        function [scaledProblem, scale] = createScaledProblem(o)
            % Create a scaled version of the problem
            scaledProblem = UtilityMaximizationProblem();
            scaledProblem.utilityFunction = o.utilityFunction;
            
            nInstruments = length(o.instruments);
            scale = ones(nInstruments,1);
            for i=1:nInstruments
                instrument = o.instruments{i};
                s = 1;
                if (isfinite(instrument.getAsk()))
                    s=instrument.getAsk();
                elseif (isfinite(instrument.getBid()) && instrument.getBid()>0)                         
                    s=instrument.getBid();
                end
                %%%%%%%%%%%%%%%%%%%%
                if 1/s ==inf
                    scale(i)=1;
                    assert(s==0);
                else
                    scale(i)=1/s;
                    assert(s>0);
                end
                
                scaledProblem.instruments{i} = RescaledInstrument( instrument, scale(i) );
                scaledProblem.payoffi{i} = scale(i)*o.payoffi{i};
            end
            
            scaledProblem.currentPosition = o.currentPosition;
            scaledProblem.scenarios = o.scenarios;
            scaledProblem.logProb = o.logProb;
            scaledProblem.payoff0 = o.payoff0;
            scaledProblem.constraints = cell( length(o.constraints));
            for i=1:length(o.constraints)
                c = o.constraints{i};
                scaledProblem.constraints{i} = c.rescale(scale);
            end
        end
        
        function [utility,quantities,qp] = optimize(o)
            % Solve the optimization problem returning the expected
            % utility and the quantities that must be held. The first
            % task is to form an equivalent problem with good scaling
            % behaviour.
%             
%                [scaledProblem, scale] = o.createScaledProblem();
%                [utility,scaledQuantities] = scaledProblem.optimizeUnscaled();
%                quantities = scaledQuantities .* scale;
             [utility,quantities,qp] = o.optimizeUnscaled();
        end
        
        function writeProblem(o, fileName)
            scaledProblem = o.createScaledProblem();
            UtilityMaximizationSolver( scaledProblem ).writeProblem(fileName);
        end
        
        function [utility,quantities,qp] = optimizeUnscaled(o)
            % Solve the optimization problem returning the expected
            % utility and the quantities that must be held. This operates
            % naively and doesn't consider rescaling
            solver = UtilityMaximizationSolverGeometric (o);
            [utility,quantities,qp] = solver.solve();
        end
        
        function utility = utilityForQuantities(o, q)
            % Compute the utility for the given quantities
            solver = UtilityMaximizationSolver(o);
            utility = solver.utilityForQuantities(q);
        end
        
        function assertConstraintsPassed(o, q, tolerance)
            % Assert all problem constraints pass within the given
            % tolerance
            solver = UtilityMaximizationSolver(o);
            solver.assertConstraintsPassed(q, tolerance);
        end
        
    end
end

