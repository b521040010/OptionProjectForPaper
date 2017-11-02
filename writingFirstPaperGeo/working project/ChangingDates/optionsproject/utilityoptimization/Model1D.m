classdef Model1D
    
    methods (Abstract)
        
        res = logPdf(model,x)
        % Compute the log pdf of the model at the point x
        
        wayPoints = getWayPoints(model)
        % Returns some standard way points for accurate numeric integration
        
        model = fit(model, S0, T, returns, frequencies)
        % Computes a model that is calibrated to the historic returns over the time
        % period. The frequencies vector is used to weight the returns so
        % that more recent data is given more weight.
    end

end
