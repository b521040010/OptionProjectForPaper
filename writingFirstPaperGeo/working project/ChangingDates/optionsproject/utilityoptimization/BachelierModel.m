classdef BachelierModel < Model1D
    % A normally distributed model for a stock for which it is possible
    % to solve the exponential utility maximization problem analytically
    
    properties
        S0;
        meanT; % The mean at time T
        sdT;   % The s.d. at time T
    end
    
    methods
        
        function model = BachelierModel()       
            % Create a model instance
            model.S0 = 1000.0;
            model.meanT = 1010.0;
            model.sdT = 200.0;
        end
        function res = pdf(o,x)
            % Compute the pdf at x
            res = exp( logPdf(o,x));
        end         
        
        function res = logPdf(o,x)
            % Compute the log of the pdf at x
            res = log(1/(o.sdT*sqrt(2*pi))) - ((x - o.meanT).^2) ./ (2* o.sdT^2);
        end          
        
        function r = mean(o)
            % Compute the mean of the model
            r = o.meanT;
        end
        
        function r = sd(o)
            % Compute the s.d. of the model
            r = o.sdT;
        end
        
        
        function wayPoints = getWayPoints(model)
        % Returns some standard way points for accurate numeric integration
            mean = model.mean();
            sd = model.sd();
            wayPoints = (-10:10)*sd + mean;
            wayPoints = wayPoints( wayPoints>0 );            
        end
                
        function u = expectedExponentialUtility( o, principal, ...
                                                 bondInvestment, ...
                                                 r, ...
                                                 T, ...
                                                 gamma ) 
            P = principal;
            lS0 = o.S0;
            lmeanT = o.meanT;
            lsdT = o.sdT;
            x = bondInvestment;
            % The expected utility in this case is calculated in the
            % Mathematica notebook exponentialBachelier. We copied the
            % result in here
            u = -((-1 + exp((gamma*(gamma*lsdT^2*(P - x)^2 - 2*lS0*(lmeanT*(P - x) + exp(r*T)*lS0*x)))/(2*lS0^2)))/gamma);
        end
        
        function bondInvestment = optimizeExponentialUtility( o, principal, ...
                                                 r, ...
                                                 T, ...
                                                 gamma ) 
            P = principal;
            lS0 = o.S0;
            lmeanT = o.meanT;
            lsdT = o.sdT;
            % See the
            % Mathematica notebook exponentialBachelier. We copied the
            % result in here
            bondInvestment = (-lmeanT*lS0 + exp(r*T)*lS0^2 + gamma*P*lsdT^2)/(gamma*lsdT^2);
        end
        
        function model = fit( bm, S0, T, returns)
            error('Not implemented');
        end
        
    end
end

