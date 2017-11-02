classdef MaxLiabilityConstraint < UtilityMaximizationConstraint
    %BoundedLiabilityConstraint Constraint that the liability must
    %   be bounded - i.e. the slope at infinity of the payoff function must
    %   be greater than or equal to 0
    
    properties (Access='private') 
        maxLoss;
    end
    
    methods 
        
        function c = MaxLiabilityConstraint( maxLoss )
            c.maxLoss = maxLoss;
        end
        
        function applyConstraint( o, utilityMaximizationSolver, separableProblem ) 
        % Apply the constraint to the separable problem
            from = utilityMaximizationSolver.nq+1;
            to = utilityMaximizationSolver.nq+utilityMaximizationSolver.nP;
            currentBound = separableProblem.blx(from:to);
            additionalBound = -o.maxLoss*ones(size(currentBound));
            separableProblem.blx(from:to) = max( currentBound, additionalBound) ;
            disp( separableProblem.blx);
        end
        
        function c = rescale( c, factor)
            % Produce a rescaled constraint when all the associated prices
            % are multiplied by the given factor            
            c.maxLoss = c.maxLoss;
        end

    end
            
end


