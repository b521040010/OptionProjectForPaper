classdef DeltaHedgeProblem
    %DELTAHEDGEPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        bsm % The associated black scholes model
        proportionTransactionCosts % The amount of transaction costs
        nPaths % the number of paths
        K % Strike of instrument hedged
        isPut % whether this is a put or a call
        r % Risk free rate
        utilityFunction
    end
    
    methods
        
        function p = DeltaHedgeProblem()
            p.bsm = BlackScholesModel();
            p.proportionTransactionCosts = 0.01;
            p.nPaths = 10000;
            p.K = p.bsm.S0;
            p.isPut = false;
            p.r = 0.10;
            p.utilityFunction = ExponentialUtilityFunction(1.0);
        end
        
        function eu = computeExpectedUtility( p, q, charge, S )
            % q is quantity sold, charge is charge for that quantity
            b = p.simulateDeltaHedge( q, charge, S );
            eu = sum(p.utilityFunction.weightedUtility(b, -log(p.nPaths) ));
        end
        
        function S= simulatePricePaths( p, nSteps ) 
            S = p.bsm.simulatePricePaths(p.nPaths, nSteps );
        end
        
        function b = simulateDeltaHedge( p, quantity, charge, S )
            % Simulate delta hedging with the given number of steps
            % returning the final bank balance
            % given the initial charge
            nSteps = size( S, 2 )-1;
            dt = p.bsm.T / nSteps;
            currS = S(:,1);
            currBsm = p.bsm;
            [~,delta] = currBsm.price( p.r, p.isPut, p.K, currS );
            b = charge - p.priceQuantity(quantity*delta,currS); % bank balance
            for i=1:nSteps-1
                currS = S(:,i+1);
                prevDelta = delta;
                currBsm.T = currBsm.T - dt;
                [~,delta] = currBsm.price( p.r,p.isPut, p.K , currS );
                b = exp(p.r * dt)*b - p.priceQuantity(quantity*(delta-prevDelta),currS);
            end
            prevDelta = delta;
            currS = S(:,nSteps+1);
            payoff = quantity*max((2*p.isPut-1)*(p.K-currS),0);
            b = exp(p.r * dt)*b - p.priceQuantity(-quantity*prevDelta,currS) - payoff;
        end
        
        function price = priceQuantity( p, quantity, askPrice )
            % Return the price for a given quantity of the stock taking
            % into account transaction costs
            bidPrice = askPrice * (1-p.proportionTransactionCosts);
            price = quantity.*((quantity>=0).*askPrice + (quantity<0).*bidPrice);
        end
        
        function price = sellersIndifferencePrice( p, quantity, priceGuess, S, passCount )
            % Perform buyers indifference pricing in the special case of
            % exponential utility and a totally liquid risk free bond
            if (nargin<5)
                passCount = 0;
            end
            V0 = 0;
            V1 = computeExpectedUtility( p, quantity, priceGuess, S );
            
            a = p.utilityFunction.gamma;
            price = exp(-p.r*p.bsm.T)*(1/a)*log((1-a*V1)/(1-a*V0)) + priceGuess;
            check = p.computeExpectedUtility( quantity, price, S );
            if (abs(check)>0.0001)
                if (passCount>3)
                    error('Unable to compute indifference price accurately');
                end
                price = p.sellersIndifferencePrice( quantity, price, S, passCount+1 );
            end
        end
        
    end
    
end

