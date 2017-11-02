function formatMathematica( m )
%Write a matrix out in mathematica format
nrows = size(m,1);
ncols = size(m,2);
fprintf('{\n');
for i=1:nrows
    fprintf('{');
    for j=1:ncols
        fprintf('%f',m(i,j));
        if j~=ncols
            fprintf(',');
        end
    end
    fprintf('}');
    if (i~=nrows)
        fprintf(',');
    end
    fprintf('\n');
end    
fprintf('}\n');
end

