classdef ArbitrageFinder2
    %ARBITRAGEFINDER2 Summary of this class goes here
    %   Detailed explanation goes here
        
    methods (Static)
        
        function [arbitrageExists, portfolio] = findArbitrage(instruments, allowBonds, allowFutures)
            portfolio = Portfolio();

            instruments = ArbitrageFinder2.filterInstruments( instruments, allowBonds, allowFutures );
            nInstruments = length(instruments);    
            
            wayPoints = unique( ArbitrageFinder2.getWayPoints(instruments) )';
            nWayPoints = length(wayPoints );
            
            payoffMatrix = zeros( nWayPoints, 2*nInstruments );
            finalDeriv = zeros(1, 2*nInstruments );
            costVector = zeros(1, 2*nInstruments );
            lb=zeros(2*nInstruments,1);
            ub=zeros(2*nInstruments,1);
            for i=1:nInstruments
                ins = instruments{i};
                payoff = ins.payoff( wayPoints )/ins.contractSize;
                ub(i)=ins.askSize;
                lb(i)=0;
                payoffMatrix(:,i) = payoff;
                payoffMatrix(:,nInstruments+i) = -payoff;
                finalDeriv(i) = ins.deltaAtInfinity();
                finalDeriv(nInstruments+i) = -ins.deltaAtInfinity();
                ub(nInstruments+i)=ins.bidSize;
                lb(nInstruments+i)=0;
                costVector(i) = ins.getAsk()/ins.contractSize;
                if isnan(costVector(i))
                    costVector(i)=Inf;
                end
                costVector(nInstruments+i) = -ins.getBid()/ins.contractSize;
                if isnan(costVector(nInstruments+i))
                    costVector(nInstruments+i)=0;
                end                
            end
            sumPayoff = sum(payoffMatrix,1);
            
            id = eye(2*nInstruments);
            
            % Constraints are:
            %  Positive quantity
            constraintMatrix1 = -id;
            constraintVector1 = zeros(2*nInstruments,1);
            %  Payoff must be positive
            constraintMatrix2 = -payoffMatrix;
            constraintVector2 = zeros( nWayPoints,1);
            %  Final deriviative must be positive
            constraintMatrix3 = -finalDeriv;
            constraintVector3 = 0;
            %  Cost must be negative
            constraintMatrix4 = costVector;
            constraintVector4 = 0;
            %  Prevent the problem being unbounded - sumPayoff<=1
            constraintMatrix5 = -sumPayoff;
            constraintVector5 = -1000;
            constraintMatrix = vertcat( constraintMatrix1, constraintMatrix2, constraintMatrix3, constraintMatrix4, constraintMatrix5);
            constraintVector = vertcat( constraintVector1, constraintVector2, constraintVector3, constraintVector4, constraintVector5);
            
            % Get the largest payoff for the least cost
            objective = sumPayoff - costVector;
            [q,~,exitFlag] = linprog( objective, constraintMatrix, constraintVector,[],[],lb,ub );            
            
            arbitrageExists = (exitFlag==1);
            if (arbitrageExists)
                errorVec = constraintMatrix * q - constraintVector;
                errorVec = errorVec(1:end-1);
                fprintf('*** Arbitrage error %d\n',max(errorVec));
                relevantInstruments = q>1e-1;
                oDash = objective(relevantInstruments);
                cMDash = constraintMatrix(:,relevantInstruments);
                cVDash = constraintVector;
                lb=lb(relevantInstruments);
                ub=ub(relevantInstruments);
                [qDash,~,exitFlag] = linprog( oDash, cMDash, cVDash,[],[],lb,ub  );     
                if (exitFlag~=1)
                    arbitrageExists = false;
                    fprintf('*** Arbitrage disappeared ***\n');
                else
                    fprintf('*** Genuine arbitrage using %d instruments ***\n', sum( relevantInstruments) );
                end
                count = 1;
                q = zeros(2*nInstruments,1);
                for i=1:(2*nInstruments)
                    if relevantInstruments(i)
                        q(i) = qDash(count);
                        count = count+1;
                    end
                end
                portfolio = Portfolio();
                qSmall = q(1:nInstruments)-q((nInstruments+1):end);
                portfolio.add( qSmall, instruments ); 
            end
            return;
            
            disp('Error');
            disp( max( constraintMatrix * q - constraintVector));
            disp('Detail');
            disp( constraintMatrix * q - constraintVector);
            
            % Two equivalent methods give different answers due to
            % approximation problems
            constraintMatrix5Dash = -ones(1,2*nInstruments);
            constraintVector5Dash = -1000000;
            objectiveDash = abs(costVector);
            constraintMatrixDash = vertcat( constraintMatrix2, constraintMatrix3, constraintMatrix4, constraintMatrix5Dash );
            constraintVectorDash = vertcat( constraintVector2, constraintVector3, constraintVector4, constraintVector5Dash);
            lb = zeros(2*nInstruments,1);            
            [~,~,exitFlagDash] = linprog( objectiveDash, constraintMatrixDash, constraintVectorDash, [], [], lb, [], q );            
%            res = msklpopt(objectiveDash,constraintMatrixDash,[],constraintVectorDash,lb,[]);
%            disp( res );
%            disp( constraintMatrixDash * q - constraintVectorDash);
%            disp( 1000000 * constraintMatrixDash * q - constraintVectorDash);
            fprintf('Exit flag vs exit flag dash (%d,%d)\n', exitFlag, exitFlagDash);
            arbitrageExists = exitFlagDash==1;
            if (arbitrageExists)
                disp('Found arbitrage error=');
                disp( max( constraintMatrixDash * q - constraintVectorDash));
                qSmall = q(1:nInstruments)-q((nInstruments+1):end);
                portfolio.add( qSmall, instruments );
                disp( portfolio.print());
            end
        end
        
        function wayPoints = getWayPoints(instruments)
            % Get the way points defined by the instruments - note that
            % we also use the way points defined by the pdf to choose the
            % integration points
            wayPoints = 0;
            nInstruments = length(instruments);
            for i=1:nInstruments
                instrument = instruments{i};
                wayPoints = horzcat(wayPoints,instrument.getWaypoints());
            end
        end
        
        function retInstruments = filterInstruments( instruments, allowBonds, allowFutures) 
            retInstruments  = cell(0,1);
            nInstruments = length( instruments );
            count = 1;
            for idx=1:nInstruments
                ins = instruments{idx};
                if isa(ins,'Bond')
                    allow = allowBonds;
                elseif isa(ins, 'Future2')
                    allow = allowFutures;
                else 
                    allow = true;
                end
                if allow
                    retInstruments{count}=ins;
                    count = count+1;    
                end
            end
        end
    end
    
end

