classdef UtilityMaximizationConstraint
    %UTILITYMAXIMIZATIONCONSTRAINT Summary of this class goes here
    %   Detailed explanation goes here
    
    methods (Abstract)        
        applyConstraint( utilityMaximizationProblem, geometricProblem )        
        % Apply the constraint to the separable problem
        
        rescale( factor)
        % Produce a rescaled constraint when all the associated prices
        % are multiplied by the given factor
    end
            
end

