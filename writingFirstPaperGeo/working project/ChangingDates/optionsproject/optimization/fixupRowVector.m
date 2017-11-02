function rowVector=fixupRowVector( name, rowVector, x)
    % Increase user friendliness by accepting row or column vectors
    % from the user supplied functions 
    % and performing error checking that length(rowVector)=length(x)
    if (size(rowVector,2)==1)
        rowVector = rowVector';
    end
    if length(rowVector)~=length(x)
        error('%s function should return same number of scenarios as it receives',name);
    end
end

