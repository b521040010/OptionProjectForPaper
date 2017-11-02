classdef QuantityConstraint < UtilityMaximizationConstraint
    %QuantityConstraint Constraint that the quantity of a given
    %   instrument must be in a certain range
    
    properties
        instrumentIndex;
        min;
        max;
    end
    
    methods 
        
        function c = QuantityConstraint( instrumentIndex, min, max)
            c.instrumentIndex = instrumentIndex;
            c.min = min;
            c.max = max;
            assert(min<=0); % Currently only supporting min<=0<=max
            assert(max>=0);
        end
        
        function applyConstraint( c, utilityMaximizationSolver, geometricProblem ) 
        % Apply the constraint to the separable problem
            s = utilityMaximizationSolver;
            
            indices = find(s.indexToInstrument==c.instrumentIndex);
            for i=indices
                short = s.indexToShort(i);
                if short
                    temp=zeros(1,s.nq+1);
                    temp(i)=1;
                    geometricProblem.addLinearUpperBound(temp,-c.min,'quantity constraint short');
                else       
                    temp=zeros(1,s.nq+1);
                    temp(i)=1;
                    geometricProblem.addLinearUpperBound(temp,c.max,'quantity constraint long');
                end
            end
            
        end
        
        function c = rescale( c, scale)
            % Produce a rescaled constraint when all the associated prices
            % are multiplied by the given factor            
            factor = scale( c.instrumentIndex );
            c.min = c.min/factor;
            c.max = c.max/factor;
        end
        
    end
            
end



