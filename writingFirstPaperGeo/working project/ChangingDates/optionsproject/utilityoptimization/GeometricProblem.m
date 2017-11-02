classdef GeometricProblem < handle
    %Mosek optimiser [res] = mskgpopt(c,a,map)
    % Min     log(w1exp(ao1.x)+w2exp(ao2.x)+...wkexp(aok.x))
    % St      aci.x <= u , i=1,2,3,...,m
    % where w=[w1 w2 ... wk]' , ao=[ao1 ; ao2 ; ... ; aok]
    % ac=[ac1 ; ac2 ; ... ; acm]'
    % lxl=size(ac,2)=nVars
    properties        
        % The number of variables
        nVars
        % inouts of the optimiser
        c
        a
        map
        blc
        buc
        % Lower bound on x variables
        blx
        % Upper bound on x variables
        bux
        
        %c=[weights;exp(1./u)]
        % coefficients of the exponential functions (objective) column
        % vector
        w
        % upper bounds of the constraints (column vector)
        u
        %a=[ao;ac]
        ao
        ac
        % Text description of constraints
         constraintDescription
        
        % Acceptable level of duality gap
        tolerance;
    end
    
    methods
        
        function o = GeometricProblem( nVars )    
            o.nVars = nVars;
            o.blx = -Inf(nVars,1);
            o.bux = Inf(nVars,1);
            %to make sure we have at least one constraint
            o.ac=zeros(1,nVars);
            o.u=zeros(1,1);
            o.constraintDescription = cell(1,0);
            o.tolerance = 1e-10;
            initMosek();     
%            o.scale=10;
        end
        
%         function writeProblem( o, filename )
%             f = fopen(filename, 'w');
%             fprintf(f,'Separable Problem\n');
%             fclose(f);
%             o.appendToFile(filename,'c');
%             dlmwrite(filename,o.c,'-append');
%             o.appendToFile(filename,'A');
%             dlmwrite(filename,o.A,'-append');
%             o.appendToFile(filename,'blc');
%             dlmwrite(filename,o.blc,'-append');
%             o.appendToFile(filename,'buc');
%             dlmwrite(filename,o.buc,'-append');
%             o.appendToFile(filename,'blx');
%             dlmwrite(filename,o.blx,'-append');
%             o.appendToFile(filename,'bux');
%             dlmwrite(filename,o.bux,'-append');
%             o.appendToFile(filename,'opr');
%             dlmwrite(filename,o.opr,'-append');
%             o.appendToFile(filename,'oprj');
%             dlmwrite(filename,o.oprj,'-append');
%             o.appendToFile(filename,'oprf');
%             dlmwrite(filename,o.oprf,'-append');
%             o.appendToFile(filename,'oprg');
%             dlmwrite(filename,o.oprg,'-append');
%         end
        
%         function appendToFile( o, filename, msg ) 
%             % Write a line of text to a file
%             f = fopen(filename, 'a');
%             fprintf(f,msg);
%             fprintf(f,'=\n');
%             fclose(f);
%         end
        
        function setObjective( o, w,ao)
%             assert( size(c,2)==1); % c should be column vector                        
%             assert( size(c,1)==o.nVars ); % Must be nVars entries in C
%             assert( sum(isfinite(c))==size(c,1))   
%             assert( sum(isfinite(oprj))==size(oprj,1))
%             assert( sum(isfinite(oprf))==size(oprf,1))
%             assert( sum(isfinite(oprg))==size(oprg,1))
            assert(size(w,2)==1);

            assert(size(ao,2)==o.nVars)
            
            assert(size(ao,1)==size(w,1))
            o.w = w;
%             o.u = u;
%             o.c = [o.w ;1./exp(o.u)];
             o.ao=ao;
%             o.ac=ac;
%             o.a = [o.ao ; o.ac];   
                     
        end
        
         function total = computeObjective( o, x )
%             % Compute the value of the objective function for a matrix
%             % of x values. Each column should be a value of x

            total=log(sum(o.w.*exp(o.ao*x),1));
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
            cValue = o.ac * x;
%             for idx = findIndices( cValue < o.u-tolerance )
%                 if (~isempty(idx))
%                     error('Lower bound on linear constraint %s failed: %d < %d\n', o.constraintDescription{idx}, cValue(idx), o.blc(idx));
%                 end
%             end
            for idx = findIndices( cValue > o.u+tolerance )
                if (~isempty(idx))
                    %error('Upper bound on linear constraint %s failed: %d > %d', o.constraintDescription{idx}, cValue(idx), o.buc(idx));
                    o.constraintDescription{idx};
                    error('Upper bound on linear constraint %s failed')
                    
                end
                
            end
        end
%         
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
%            assert( sum(isfinite(b))==size(b,1));% No NaN values
            
            currentNConstraints = size(o.ac,1);
             
            if (upperBound)
                 o.ac = vertcat( o.ac, A );
                 o.u = vertcat( o.u, b );
            end
            if (lowerBound)
                 o.ac = vertcat( o.ac, -A );
                 o.u = vertcat( o.u, -b );
            end
            nConstraints = size(A,1);
            for i=1:nConstraints
                o.constraintDescription{currentNConstraints+i} = sprintf('%s - row %d',desc,i);
            end
        end
        
        
        function [objective, x, res] = optimize(o) 
            for i=1:length(o.u)
                if(o.u(i)>0)
                    o.ac(i,:)=o.ac(i,:)./o.u(i);
                    o.u(i)=1;
                end
            end
            
%             o.c = [o.w ;1./exp(o.u/o.scale)];
%             o.a = [o.ao ; o.ac/o.scale];   
            o.c = [o.w ;1./exp(o.u)];
            o.ac=sparse(o.ac);
            o.a = [o.ao ; o.ac]; 
            size(o.a);
            temp=1:1:size(o.ac,1);
            o.map = [zeros(size(o.w,1),1);temp'];
            assert(size(o.map,2)==1); 
            assert(size(o.u,2)==1);
            assert(size(o.ac,2)==o.nVars);
%             if (size(o.A,2))==0
%                 error('You must add at least one linear constraint. This is because Mosek uses the constraint matrix A to determine the problem size');
%             end
%             
            param = [];            
           %  param.MSK_DPAR_INTPNT_NL_TOL_REL_GAP = 10^(-16);
            % param.MSK_DPAR_INTPNT_CO_TOL_DFEAS=10^(-16);
           % param.MSK_IPAR_INTPNT_MAX_ITERATIONS = 100000000;
            param.MSK_IPAR_LOG = 0;
%            toc
            tic
            
            [res] = mskgpopt(o.c,o.a,o.map); 
            toc
%             if (~strcmp(res.sol.itr.solsta,'OPTIMAL') ...
%                 && ~strcmp(res.sol.itr.solsta,'NEAR_OPTIMAL'))
%                 % Repeat the optimization with debug information
%                 param.MSK_IPAR_LOG = 10;
%                 mskscopt(o.opr,o.opri,o.oprj,o.oprf,o.oprg,o.c,o.A,o.blc,o.buc,o.blx, o.bux, param);                                     
%                 error('Unable to find solution to optimization problem: %s',res.sol.itr.solsta);
%             end
%             x = res.sol.itr.xx;
%             objective = res.sol.itr.pobjval;
     
            x=res.sol.itr.xx;
            objective = o.computeObjective(x);
        end
    end
    
%     methods (Static)
%         function term = sepEval( oprName, xj, f, g )
%         % Evaluate a separable function
%             switch oprName
%                 case 'ent'
%                     term = f .* xj .* log(xj );
%                 case 'exp'
%                     term = f .* exp( g.* xj );
%                 case 'log'
%                     term = f .* log( xj );
%                 case 'pow'
%                     term = f .* xj.^g;
%                 otherwise
%                     error('Invalid function: %s',oprName);
%             end
%         end            
%     end
%     
end

function ret = findIndices( vec )
indices = find(vec, 1);
ret = reshape(indices, [1 length(indices )]);
end