classdef StudentTModel < Model1D
    % Model based on the future index prices being tdistributed
    % distributed with center parameter mu, scale parameter sigma
    % and degrees of freedom parameter nu
    
    properties
        S0;
        nu;  % degrees of freedom
        mu;   % centering parameter
        sigma; % sigma parameter
        T;
    end
    
    methods
        function model = StudentTModel()     
            model.S0=1000;
            model.nu = 3;
            model.mu = 0;
            model.sigma = 0.2;
            model.T = 1/12;
        end
       
        
        function res = pdf(m,x)
        % Compute the pdf of the model
            s = log(x);
            res = gamma((m.nu+1)/2)/(gamma(m.nu/2)*sqrt(pi*m.nu)*m.sigma) ...
                *( 1 + 1/m.nu*((s-m.mu)/m.sigma).^2).^(-(m.nu+1)/2) ...
                ./ x;
        end         
        
        function res = logPdf(m,x)
        % Compute the log pdf of the model
            s = log(x);
            res = log(gamma((m.nu+1)/2)/(gamma(m.nu/2)*sqrt(pi*m.nu)*m.sigma)) ...
                + (-(m.nu+1)/2).*log( 1 + 1/m.nu*((s-m.mu)/m.sigma).^2) ...
                - s;
        end          
        
        function [ price] = price(o, r, isPut, strike, spot)
            %price call option
            if isPut==1
                error('the calibration for a put has not been implimented')
            end
            
            assert( islogical( isPut ));
            vol=o.sigma;
            nuu=o.nu;
            TT=o.T;
            if (nargin<5)
                spot = o.S0;
            end    
    
            function ret=integrand(x)
%                s = log(x);
                pdf = gamma((nuu+1)/2)/(gamma(nuu/2)*sqrt(pi*nuu)*vol) ...
                        *( 1 + 1/nuu*((log(x)-(log(spot)+r))/vol).^2).^(-(nuu+1)/2) ...
                        ./ x;
%                 pdf = gamma((nuu+1)/2)./(gamma(nuu/2).*x*sqrt(pi*nuu)*vol).*( 1 + 1/nuu*((s-(log(spot)+r))/vol).^2).^(-(nuu+1)/2) 
       
%                payoff=((1-isPut)*max(x-strike,0)+isPut*max(strike-x,0));
               payoff=x-strike;
%                 ret=payoff.*exp(-r*TT).*o.pdf(x);.
                ret=payoff.*exp(-r*TT).*pdf;
            end
            price=integral(@integrand,strike,100000,'RelTol',0,'AbsTol',1e-12);
        end

            function vol = impliedSigma( m, r, strike, isPut, optionPrice)
                %compute implied sigma and nu when the index is log student t
                %distributed
                %we use 2 in-the-money call options to calibrate the parameters
                %optionPrice=optionPrice/100;
                function F = root2d(x)
                    %x(1)=sigma, x(2)=nu
                    m.sigma=x;
                    F=m.price( r, isPut, strike )-optionPrice;
                end
                fun = @root2d;
                options = optimoptions('fsolve','MaxFunctionEvaluations',1000,'OptimalityTolerance',10^(-10));
                [vol,~, ret] = fsolve(fun,0.1,options);
                if vol>1 || vol<=0
                    error('calibrated sigma is strange')
                end
                m.price( r, isPut, 2440 )
%                 if ret<=0
%                     fprintf('Cannot find implied volatility\n');
%                     fprintf( 'r=%f, T=%f, S0=%f, strike=%f, isPut=%f, optionPrice%f\n', r, m.T, m.S0, strike, isPut, optionPrice );
%                     error('Failed');
%                 end         

            end            
        
        
            function nu = impliedNu( m, r, strike, isPut, optionPrice,vol)
                %compute implied sigma and nu when the index is log student t
                %distributed
                %we use 2 in-the-money call options to calibrate the parameters
                %optionPrice=optionPrice/100;

                function F = root2d(x)
                    %x(1)=sigma, x(2)=nu
                    m.nu=x;
                    F=m.price( r, isPut, strike )-optionPrice;
                end
                fun = @root2d;
                options = optimoptions('fsolve','MaxFunctionEvaluations',1000,'OptimalityTolerance',10^(-10),'Display', 'off');
                [nu,~, ret] = fsolve(fun,5,options);
                if ret<=0
                    fprintf('Cannot find implied volatility\n');
                    fprintf( 'r=%f, T=%f, S0=%f, strike=%f, isPut=%f, optionPrice%f\n', r, m.T, m.S0, strike, isPut, optionPrice );
                    error('Failed');
                end         
                if nu<=0
                    error('calibrated nu is negative')
                end
            end

            function sigma = impliedVolatilityBS( m, r, strike, isPut, optionPrice )
                % Compute the implied volatility of an option
                %optionPrice=optionPrice/100;
                b = BlackScholesModel();
                b.T=m.T;
                b.S0=m.S0;
                function priceDiff = f( sigma )
                    b.sigma = sigma;
                    priceDiff = b.price( r, isPut, strike ) - optionPrice;
                end
                options = optimset('fsolve');
                options = optimset(options, 'Display', 'off');
                [sigma,~, ret] = fsolve( @f, 0.3, options );
                if ret<=0
                    fprintf('Cannot find implied volatility\n');
                    fprintf( 'r=%f, T=%f, S0=%f, strike=%f, isPut=%f, optionPrice%f\n', r, m.T, m.S0, strike, isPut, optionPrice );
                    error('Failed');
                end
            end            
            
        function m = fit( m, S0, T, returns, weights)
            % Fit the model to historic return data            
            fittedDist = fitdist( log(returns+1), 'tlocationscale', 'frequency', weights);
            m.mu = log(S0)+fittedDist.mu;
            m.sigma = fittedDist.sigma;
            m.nu = fittedDist.nu;
            m.T = T;

        end

        function [S, times] = simulatePricePaths( model, nPaths, nSteps )
            if nSteps>1
                error('the code for nSteps > 1 has not been written.')
            end
            %dt = model.T/nSteps;
             p = haltonset(1,'Skip',2);
            rn=net(p,nPaths);
            s=icdf('tLocationScale',rn,model.mu,model.sigma,model.nu);
            
            S = horzcat(model.S0*ones(nPaths,1),exp(s));
            times = linspace(0,model.T,nSteps+1);
        end
        
        function returns=getReturns(m)
            %open file SPX Index and get return form (log(ST/S0))
            index=xlsread('../SPX Index');
            n=length(index);
            returns=(index(2:n)-index(1:n-1))./index(1:n-1);
        end
        
        function wayPoints = getWayPoints(m)
        % Returns some standard way points for accurate numeric integration
            beforeScaling = tinv(0.01:0.01:0.99,m.nu);
            wayPoints = [0 exp(beforeScaling * m.sigma + m.mu)];
        end
        
    end
end

