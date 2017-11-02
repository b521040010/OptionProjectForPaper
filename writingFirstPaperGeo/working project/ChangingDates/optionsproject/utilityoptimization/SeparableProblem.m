classdef SeparableProblem < handle
    %SEPARABLEPROBLEM Represents a separable optimization problem
    %  It is intended to make mskscopt a little more practical to use
    
    properties        
        % The number of variables
        nVars
        % Linear part of objective definition
        c
        % Linear part of constraint definition
        A
        blc
        buc
        % Lower bound on x variables
        blx
        % Upper bound on x variables
        bux
        % Definition of non linear functions
        opr
        opri
        oprj
        oprf
        oprg
        
        % Text description of constraints
        constraintDescription
        
        % Acceptable level of duality gap
        tolerance;
    end
    
    methods
        
        function o = SeparableProblem( nVars )    
            o.nVars = nVars;
            o.blx = -Inf(nVars,1);
            o.bux = Inf(nVars,1);
            o.constraintDescription = cell(1,0);
            o.tolerance = 1e-10;
            initMosek();            
        end
        
        function writeProblem( o, filename )
            f = fopen(filename, 'w');
            fprintf(f,'Separable Problem\n');
            fclose(f);
            o.appendToFile(filename,'c');
            dlmwrite(filename,o.c,'-append');
            o.appendToFile(filename,'A');
            dlmwrite(filename,o.A,'-append');
            o.appendToFile(filename,'blc');
            dlmwrite(filename,o.blc,'-append');
            o.appendToFile(filename,'buc');
            dlmwrite(filename,o.buc,'-append');
            o.appendToFile(filename,'blx');
            dlmwrite(filename,o.blx,'-append');
            o.appendToFile(filename,'bux');
            dlmwrite(filename,o.bux,'-append');
            o.appendToFile(filename,'opr');
            dlmwrite(filename,o.opr,'-append');
            o.appendToFile(filename,'oprj');
            dlmwrite(filename,o.oprj,'-append');
            o.appendToFile(filename,'oprf');
            dlmwrite(filename,o.oprf,'-append');
            o.appendToFile(filename,'oprg');
            dlmwrite(filename,o.oprg,'-append');
        end
        
        function appendToFile( o, filename, msg ) 
            % Write a line of text to a file
            f = fopen(filename, 'a');
            fprintf(f,msg);
            fprintf(f,'=\n');
            fclose(f);
        end
        
        function setObjective( o, c, opr, oprj, oprf, oprg )
            assert( size(c,2)==1); % c should be column vector                        
            assert( size(c,1)==o.nVars ); % Must be nVars entries in C
            assert( sum(isfinite(c))==size(c,1))   
            assert( sum(isfinite(oprj))==size(oprj,1))
            assert( sum(isfinite(oprf))==size(oprf,1))
            assert( sum(isfinite(oprg))==size(oprg,1))
            
            o.c = c;
            o.opr = opr;                        
            % Validate the nonlinear function
            assert( size(opr, 2)==3 ); % opr should contain operation names
            n = size(opr,1);
            assert( size(oprj, 2)==1 ); 
            assert( size(oprf, 2)==1 );            
            assert( size(oprg, 2)==1 );
            assert( size(oprj, 1)==n );
            assert( size(oprf, 1)==n );
            assert( size(oprg, 1)==n );            
            o.opri = zeros(n,1);
            o.oprj = oprj;
            o.oprf = oprf;
            o.oprg = oprg;
        end
        
        function total = computeObjective( o, x )
            % Compute the value of the objective function for a matrix
            % of x values. Each row should be a value of x
            n = size(o.opr,1);
            total = o.c' * x';
            for idx=1:n
                oprName = o.opr(idx,:);
                j = o.oprj(idx,1);
                xj = x(:,j);
                f = o.oprf(idx,1);
                g = o.oprg(idx,1);
                term = SeparableProblem.sepEval( oprName, xj, f, g );
                total = total + term;
            end            
        end
        
        function assertConstraintsPassed( o, x, tolerance )
            % Confirm that the constraints are passwed within a given 
            % tolerance
            for idx = findIndices( x < o.blx-tolerance )
                if (~isempty(idx))
                    error('Lower bound constraint %d failed: %d < %d\n', idx, x(idx), o.blx(idx));
                end
            end
            for idx = findIndices( x > o.bux+tolerance )
                disp('bux');
                disp(o.bux);
                if (~isempty(idx))
                    error('Upper bound constraint %d failed: %d > %d\n', idx, x(idx), o.bux(idx));
                end
            end
            cValue = o.A * x;
            for idx = findIndices( cValue < o.blc-tolerance )
                if (~isempty(idx))
                    error('Lower bound on linear constraint %s failed: %d < %d\n', o.constraintDescription{idx}, cValue(idx), o.blc(idx));
                end
            end
            for idx = findIndices( cValue > o.buc+tolerance )
                if (~isempty(idx))
                    error('Upper bound on linear constraint %s failed: %d > %d', o.constraintDescription{idx}, cValue(idx), o.buc(idx));
                end
            end
        end
        
        function addLinearConstraint( o, A, b, desc ) 
            % Add a constraint of the form Ax == b;
            addLinearConstraintPrivate(o,A,b,desc,true,true);
        end

        function addLinearUpperBound( o, A, b, desc ) 
            % Add a constraint of the form Ax < b;
            addLinearConstraintPrivate(o,A,b,desc,false,true);
        end

        function addLinearLowerBound( o, A, b, desc ) 
            % Add a constraint of the form Ax > b;
            addLinearConstraintPrivate(o,A,b,desc,true,false);
        end

        
        function addLinearConstraintPrivate( o, A, b, desc, lowerBound, upperBound )
            % Add a constraint of the form Ax = b;
            assert( size(A,2)==o.nVars); % number of variables should match
            assert( size(b,2)==1); % b should be a column vector            
            assert( size(A,1)==size(b,1)); % number of variables should match
            
            assert( sum(sum( isfinite(A)))==size(A,1)*size(A,2)); % No NaN values
            assert( sum(isfinite(b))==size(b,1));% No NaN values
            
            currentNConstraints = size(o.A,1);
            o.A = vertcat( o.A, A );
            if (lowerBound)
                o.blc = vertcat( o.blc, b );
            else
                o.blc = vertcat( o.blc, -Inf( size(b)));
            end
            if (upperBound)
                o.buc = vertcat( o.buc, b );            
            else
                o.buc = vertcat( o.buc, Inf( size(b)));
            end
            nConstraints = size(A,1);
            for i=1:nConstraints
                o.constraintDescription{currentNConstraints+i} = sprintf('%s - row %d',desc,i);
            end
        end
        
        
        function [objective, x, res] = optimize(o) 
            if (size(o.A,2))==0
                error('You must add at least one linear constraint. This is because Mosek uses the constraint matrix A to determine the problem size');
            end
            
            param = [];            
             param.MSK_DPAR_INTPNT_NL_TOL_REL_GAP = 10^(-3);
             param.MSK_DPAR_INTPNT_CO_TOL_DFEAS=10^(-12);
            param.MSK_IPAR_INTPNT_MAX_ITERATIONS = 100000000;
            param.MSK_IPAR_LOG = 0;
            [res] = mskscopt(o.opr,o.opri,o.oprj,o.oprf,o.oprg,o.c,o.A,o.blc,o.buc,o.blx, o.bux, param);                         
%             if (~strcmp(res.sol.itr.solsta,'OPTIMAL') ...
%                 && ~strcmp(res.sol.itr.solsta,'NEAR_OPTIMAL'))
%                 % Repeat the optimization with debug information
%                 param.MSK_IPAR_LOG = 10;
%                 mskscopt(o.opr,o.opri,o.oprj,o.oprf,o.oprg,o.c,o.A,o.blc,o.buc,o.blx, o.bux, param);                                     
%                 error('Unable to find solution to optimization problem: %s',res.sol.itr.solsta);
%             end
            x = res.sol.itr.xx;
            objective = res.sol.itr.pobjval;
        end
    end
    
    methods (Static)
        function term = sepEval( oprName, xj, f, g )
        % Evaluate a separable function
            switch oprName
                case 'ent'
                    term = f .* xj .* log(xj );
                case 'exp'
                    term = f .* exp( g.* xj );
                case 'log'
                    term = f .* log( xj );
                case 'pow'
                    term = f .* xj.^g;
                otherwise
                    error('Invalid function: %s',oprName);
            end
        end            
    end
    
end

function ret = findIndices( vec )
indices = find(vec, 1);
ret = reshape(indices, [1 length(indices )]);
end