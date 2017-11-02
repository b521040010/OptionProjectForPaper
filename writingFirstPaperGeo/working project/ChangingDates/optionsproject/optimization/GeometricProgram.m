classdef GeometricProgram
    
    properties
        a;
        c;
        map;
    end
    
    methods                 
        
        function gp = setObjective(gp, aObj, cObj)
            % Set the objective function
            gp.a = aObj;
            gp.c = cObj;
            gp.map = zeros(size( aObj,1 ),1);
        end
        
        
        function gp = addConstraint(gp, aCons, cCons)
            % Add a single constraint
            nrows = size(aCons,1);
            gp = gp.addConstraints( aCons, cCons, ones( nrows, 1 ) );
        end
        
        function gp = addConstraints(gp, aCons, cCons, map)
            % Add multiple constraints
            gp.a = [gp.a ; aCons];
            gp.c = [gp.c ; cCons];
            maxMap = max( gp.map );
            assert( maxMap>=0);
            gp.map = [gp.map ; map+maxMap];
        end
        
    end


end

