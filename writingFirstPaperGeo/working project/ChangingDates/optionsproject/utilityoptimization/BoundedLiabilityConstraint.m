classdef BoundedLiabilityConstraint < UtilityMaximizationConstraint
    %BoundedLiabilityConstraint Constraint that the liability must
    %   be bounded - i.e. the slope at infinity of the payoff function must
    %   be greater than or equal to 0
    
    methods 
        
        function b = BoundedLiabilityConstraint()
        end
        
        function applyConstraint( o, utilityMaximizationSolver, separableProblem ) 
        % Apply the constraint to the separable problem
            nInstruments = utilityMaximizationSolver.nq;
            vec = zeros( 1, utilityMaximizationSolver.nq +  utilityMaximizationSolver.nP );
            for i=1:nInstruments
                iIdx = utilityMaximizationSolver.indexToInstrument(i);
                short = utilityMaximizationSolver.indexToShort(i);
                instrument = utilityMaximizationSolver.p.instruments{iIdx};
                delta = instrument.deltaAtInfinity();
                if (short)
                    vec(i)= delta;
                else
                    vec(i)=-delta;
                end
            end
            separableProblem.addLinearLowerBound(vec,0, 'Bounded liability constraint');
        end
        
        function c = rescale( c, factor)
            % Produce a rescaled constraint when all the associated prices
            % are multiplied by the given factor            
        end

    end
            
end

