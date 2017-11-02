function [ root,quantities ] = lineSearch( f, estimate, tolerance )
%LINESEARCH Find a root of an increasing function
%this line search was changed for portfolio optimization. The output
%includes quantities. N.B. to have quantities, we also change the function
%f in indifferenceprices


initialWidth = abs(estimate)/100.0;
if (initialWidth<1e-6) 
    initialWidth = 1e-6;
end

width = initialWidth;
lowerVal = Inf;
cont = 1;
while (cont)
    lowerGuess = estimate - width;
    [newLowerVal,~] = f(lowerGuess);
%    if newLowerVal>lowerVal
%        error('The function f is non increasing');
%    end
    lowerVal = newLowerVal;
    if (lowerVal<tolerance)
        cont = false;
    end
    width = width*2;
end

width = initialWidth;
upperVal = -Inf;
cont = 1;
while (cont)
    upperGuess = estimate + width;
    [newUpperVal,~] = f(upperGuess);
%    if newUpperVal<upperVal
%        error('The function f is non increasing');
%    end
    upperVal = newUpperVal;
    if (upperVal>tolerance)
        cont = false;
    end
    width = width*2;
end

[root,quantities] = solveByBisection( f, lowerGuess, upperGuess, lowerVal, upperVal, tolerance );

end

function [r,quantities] = solveByBisection( f, lowerGuess, upperGuess, lVal, uVal, tolerance )
% Find a root of an increasing functioni by bisection
middle = 0.5*(lowerGuess+upperGuess);
lowerGuess
upperGuess
[utility,quantities] = f(middle);
mVal=utility
tolerance=0.000001
%if ( abs(mVal)<tolerance  && abs(lowerGuess-upperGuess)<tolerance )
if ( abs(mVal)<0.1  && abs(lowerGuess-upperGuess)<tolerance )
    r = middle;
elseif mVal>0
    r = solveByBisection(f,lowerGuess,middle,lVal,mVal,tolerance);
else
    r = solveByBisection(f,middle,upperGuess,mVal,uVal,tolerance);
end    
end