classdef UtilityMaximizationProblem1D < matlab.mixin.Copyable
    % UtilityMaximizationProblem1D  a 1 dimensional utility maximization
    %      problem that uses a 1-d Guassian quadrature approach to compute
    %      integrals.
    %
    
    properties %(Access='private')               
        
        instruments % The instruments available, stored as cell array
        % The current position represented as an instrument object
        currentPosition
        
        model % A 1-d model
                
        % The utility function
        utilityFunction
        
        % Additional constraints
        constraints;
        
        % Object we delegate to
        delegate;
        
        % tolerance for duality gap
        tolerance;
    end
    
    methods
        
        function o = UtilityMaximizationProblem1D()
            o.utilityFunction = PowerUtilityFunction(1);
            o.currentPosition = Bond(0,0,0,0,Inf,Inf);
            o.delegate = [];
            o.tolerance = 1e-9;
        end
        
        function setTolerance( o, tol )
            o.tolerance = tol;
        end
        
        function instruments = getInstruments(o)
            instruments = o.instruments;
        end
        
         function fig = plotPortfolio( o, t, quantities, initialWealth )
            portfolio = Portfolio();
            portfolio.add( quantities, o.getInstruments() );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            totalInvest=0;
            for i = 1:length(portfolio.map)
                if quantities(i)>=0
                    totalInvest = totalInvest+quantities(i)*o.instruments{i}.getAsk();
                else
                    totalInvest = totalInvest+quantities(i)*o.instruments{i}.getBid();
                end
            end
            totalInvest
            temp=values(o.currentPosition.map);
%             initialWealth=temp{1}.quantity;
            instrument=temp{1}.instrument;
%             temp=portfolio.map.values
%             temp{2}
            %We put money we have left from investing to the bank account
            %**initialWealth-totalInvest**
            %We substract the total value of the portfolio by initialWealth
            %so that we can get net profit at T
            %portfolio.add([initialWealth-totalInvest-initialWealth],{instrument});
%            assert(initialWealth-totalInvest>=0);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             cellI{1}=o.currentPosition;
%             portfolio.add( 1, cellI );
            wayPoints = o.getWayPoints();
            wayPoints = sort( unique(wayPoints) );
            range = wayPoints(end)-wayPoints(1);
            startPoint = wayPoints(1)-range/20;
            endPoint = wayPoints(end)+range/20;
            wayPoints = [ startPoint wayPoints endPoint ];
%            payoffs = portfolio.payoff( wayPoints );
            payoffs = portfolio.payoff( wayPoints )+initialWealth;
%            fig = figure();
            hold on;
            plot( wayPoints, payoffs);
%             hold on;
%             plot( wayPoints,0);

            title(t);
            ylabel('Payoff');
            xlabel('Stock price');
        end
        
        function [p,quantities] = indifferencePrice( o, cashInstrument, quantity, instrument, guess )
            % Compute the indifference price of the given quantity
            % of the instrument
            [p,quantities] = o.createDelegate().indifferencePrice(cashInstrument, quantity, instrument, guess );
        end
        
        function p = indifferencePriceExponentialUtility( o, cashInstrument, quantity, instrument, guess )
            % Compute the indifference price of the given quantity
            % of the instrument
            p = o.createDelegate().indifferencePriceExponentialUtility(cashInstrument, quantity, instrument, guess);
        end
        
        
        function addConstraint( o, c )
            idx = length( o.constraints )+1;
            o.constraints{idx} = c;
            o.delegate = [];
        end
        
        function setModel(o, model )
            o.model = model;
            o.delegate = [];
        end
        
        function model = getModel(o)
            model = o.model;
        end
        
        function setCurrentPosition( o, currentPosition ) 
            % Sets the current position
            o.currentPosition = currentPosition;
            o.delegate = [];
        end
        
        function addInstrument( o, i )
            % Adds a possible instrument
            idx = length( o.instruments )+1;
            o.instruments{ idx} = i;
            o.delegate = [];
        end
        
        function removeInstrument( o, instrument )
            instrumentText = instrument.print();
            function accept = filter( otherInstrument ) 
                otherInstrumentText = otherInstrument.print();
                accept = ~strcmp(instrumentText,otherInstrumentText);                
            end
            removedCount = o.filterInstruments( @filter );
            assert( removedCount==1);
        end
        
        function removedCount = filterInstruments( o, filter )
            % Filters the instruments only keeping those passing the filter
            n = length( o.instruments );
            mask = zeros(n,1);
            for i=1:n
                mask(i) = filter( o.instruments{i} );                
            end
            newN = sum( mask );
            oldInstruments = o.instruments;
            o.instruments = cell(1,newN);
            countL = 1;
            countR = 1;
            for i=1:n
                if mask(i)
                    o.instruments{countL} = oldInstruments{countR};
                    countL = countL+1;
                end
                countR = countR+1;
            end
            o.delegate = [];
            removedCount = length( oldInstruments ) - newN;
        end
        
        function setUtilityFunction( o, utilityFunction )
            o.utilityFunction = utilityFunction;
            o.delegate = [];
        end        
        
        function ret = getUtilityFunction(o)
            ret = o.utilityFunction;
        end
            

        function [utility,quantities,qp] = optimize(o)
            % Solve the optimization problem returning the expected
            % utility and the quantities that must be held. The first
            % task is create a multi dimensional optimizer and set its
            % quad rule
            o.delegate = createDelegate(o);
            [utility,quantities,qp ] = o.delegate.optimize();
        end
        
        
        function utility = utilityForQuantities(o, q)
            % Assert all problem constraints pass within the given
            % tolerance
            o.delegate = createDelegate(o);
            utility = o.delegate.utilityForQuantities(q);
        end
        
        function assertConstraintsPassed(o, q, tolerance)
            % Assert all problem constraints pass within the given
            % tolerance
            o.delegate = createDelegate(o);
            o.delegate.assertConstraintsPassed(q, tolerance);
        end

        function wayPoints = getWayPoints(o)
            % Get the way points defined by the instruments - note that
            % we also use the way points defined by the pdf to choose the
            % integration points
            wayPoints = [];
            nInstruments = length(o.instruments);
            for i=1:nInstruments
                instrument = o.instruments{i};
                wayPoints = horzcat(wayPoints,instrument.getWaypoints());
            end
            wayPoints = horzcat( wayPoints, o.currentPosition.getWaypoints());            
        end
        
        function writeProblem(o,fileName)
            d = o.createDelegate();
            d.writeProblem(fileName);
        end
    
        function delegate = createDelegate(o)
            if (~isempty(o.delegate))
                delegate = o.delegate;
            else

                % Create a multi dimensional optimizer
                delegate= UtilityMaximizationProblem();

                %Choose a good quadrature rule
                wayPoints = o.getWayPoints();
                wayPoints = 5:5:5000;
                quadRule = QuadRule.adapted( @(x) o.model.pdf(x), wayPoints, o.model.getWayPoints() );
                x = quadRule.x;
                lp = log(quadRule.weights) + o.model.logPdf(x);
                j=1;
                for i=1:length(x)
                if exp(lp(i))>0
                    xx(j)=x(i);
                    lpp(j)=lp(i);
                    j=j+1;
                end
                end
                delegate.setQuadRule( xx', lpp' );  
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %Monte Carlo
                %For importance sampling
%                 if strcmp(class(o.model),'BlackScholesModel')
%                     disp('This is BlackScholesModel');
%                 modelT = StudentTModel();
%                 modelT.S0=o.model.S0;
%                 modelT.nu = 5; 
%                 mu = (o.model.mu - 0.5*o.model.sigma.^2) * o.model.T + log(o.model.S0);
%                 sigma = o.model.sigma*sqrt(o.model.T);       
%                 modelT.mu = mu;            
%                 modelT.sigma = sigma;
%                 modelT.T = o.model.T;
%                 prices = modelT.simulatePricePaths(50000,1);
%                 scenarios = prices(:,end); 
%                 scenarios = sort(scenarios);
%                 wayPoints = scenarios;
%                 x = wayPoints;
%                 weights = (1/length(wayPoints))*ones(1,length(wayPoints));
%                 lp = log(weights)'+o.model.logPdfImpSamp(x);
%                 delegate.setQuadRule( x, lp );
%                 %%%%%%%%%%%%
%                 
%                 else
%                     disp('This is not BlackScholesModel')
%                 prices = o.model.simulatePricePaths(50000,1);
%                 scenarios = prices(:,end); 
%                 scenarios = sort(scenarios);
%                 wayPoints = scenarios;
%                 x = wayPoints;
%                 weights = (1/length(wayPoints))*ones(1,length(wayPoints));
%                 lp = log(weights)';
%                 delegate.setQuadRule( x, lp );
%                 end
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                

                % Everything else is straightforward
                nInstruments = length(o.instruments);                
                for i=1:nInstruments
                    delegate.addInstrument(o.instruments{i});
                end
                delegate.setCurrentPosition( o.currentPosition );
                delegate.setUtilityFunction( o.utilityFunction );
                nConstraints = length(o.constraints);
                for i=1:nConstraints
                    delegate.addConstraint(o.constraints{i});
                end
                delegate.tolerance = o.tolerance;
            end        
        end
        
        function [sup,sub] = superHedgePrice( o, claim ,quantity )
            %all options available must be used since the waypoints are
            %used as the prices that the underlying can go
            p = Portfolio();
            p.add([quantity],{claim})
            sup = o.superHedgePricePrivate( p);
            p = Portfolio();
            p.add([-quantity],{claim})
            sub = -o.superHedgePricePrivate( p );
%             p = Portfolio();
%             q(1) = -quantity;
%             i{1} = claim;
%             p.add( q, i );
%             sub = -o.superHedgePricePrivate( p );
        end
        
        function [sup] = superHedgePricePrivate( o, claim )
            % Compute the super and sub hedging price for the instrument
            wayPoints = unique(horzcat(0, o.getWayPoints(), claim.getWaypoints() ));
            %wayPoints = [wayPoints];
            wayPoints = [5:5:6000];
            nWayPoints = length( wayPoints );
            
            nInstruments = length( o.instruments );
            costVector = zeros( 1, nInstruments*2+1);
            wayPointPayoffs = zeros( nWayPoints, nInstruments*2+1);
            finalGradients = zeros( 1, nInstruments*2+1);
            
            objective = zeros( 1, nInstruments*2+1);
            objective(2*nInstruments+1)=1;
            costVector(2*nInstruments+1)=-1;
            
            ub=zeros(2*nInstruments+1,1);
            ub(2*nInstruments+1,1)=Inf;
            for i=1:nInstruments
                instrument = o.instruments{i};
                pp = instrument.payoff( wayPoints );
                costVector(i)=instrument.getAsk();
                costVector(i+nInstruments)=-instrument.getBid();
                wayPointPayoffs(:,i)=pp';
                wayPointPayoffs(:,i+nInstruments)=-pp';
                finalGradients(i)=instrument.deltaAtInfinity();
                finalGradients(nInstruments+i)=-instrument.deltaAtInfinity();
                ub(i,1)=instrument.getAskSize;
                ub(i+nInstruments,1)=instrument.getBidSize;

            end
            
                        
            minPayoff = claim.payoff( wayPoints )';
            minDeltaAtInfinity = claim.deltaAtInfinity();
            
            % Positive quantity constraint
            positivityMatrix = -horzcat(eye(2*nInstruments), zeros(2*nInstruments,1));
            positivityVector = zeros(2*nInstruments,1);
            
            constraintMatrix = vertcat( positivityMatrix, costVector, -wayPointPayoffs, -finalGradients );
            constraintVector = vertcat( positivityVector, 0, -minPayoff, -minDeltaAtInfinity );
%             constraintMatrix = vertcat( positivityMatrix, costVector, -wayPointPayoffs  );
%             constraintVector = vertcat( positivityVector, 0, -minPayoff);    

            o.instruments;
            %objective
            %constraintMatrix
            %constraintVector
            [x, sup, exitFlag] = linprog( objective, constraintMatrix, constraintVector,[],[],[],ub );

            if (exitFlag==-3)
                minPayoff = zeros(size(minPayoff));
                minDeltaAtInfinity = 0;
                constraintMatrix1 = vertcat( positivityMatrix, costVector, -wayPointPayoffs, -finalGradients, objective, -objective );
                constraintVector1 = vertcat( positivityVector, 0, -minPayoff, -minDeltaAtInfinity, 2e6, 2e6  );
                [q1, ~, ~] = linprog( objective, constraintMatrix1, constraintVector1 );            
                disp('wayPoints');
                disp(wayPoints);
                disp('Arbitrage portfolio');
                p = Portfolio();
                q = q1(1:nInstruments)-q1(nInstruments+1:2*nInstruments);
                p.add(q,o.instruments);
                disp(p.print());
                disp('Cost of portfolio');
                disp( costVector(1:(end-1))*q1(1:(end-1)) );
                disp('Payoffs + final gradient');
                disp((wayPointPayoffs*q1(1:end))');
                error('Arbitrage found');
            end
            if ~(exitFlag==1 || exitFlag==-2)
                error('Unexpected exit flag from linprog %d\n', exitFlag );
            end
            % If not feasible then infinite super hedging costs
            if (exitFlag==-2) 
                sup = Inf;
            end
            
        end
        
        
    end
end

