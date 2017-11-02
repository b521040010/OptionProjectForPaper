function [ root ] = lineSearch( f, estimate, tolerance )
%LINESEARCH Find a root of an increasing function


initialWidth = abs(estimate)/100.0;
if (initialWidth<1e-6) 
    initialWidth = 1e-6;
end

width = initialWidth;
lowerVal = Inf;
cont = 1;
while (cont)
    lowerGuess = estimate - width;
    newLowerVal = f(lowerGuess);
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
    newUpperVal = f(upperGuess);
%    if newUpperVal<upperVal
%        error('The function f is non increasing');
%    end
    upperVal = newUpperVal;
    if (upperVal>tolerance)
        cont = false;
    end
    width = width*2;
end

root = solveByBisection( f, lowerGuess, upperGuess, lowerVal, upperVal, tolerance );

end

function r = solveByBisection( f, lowerGuess, upperGuess, lVal, uVal, tolerance )
% Find a root of an increasing functioni by bisection
middle = 0.5*(lowerGuess+upperGuess);
mVal = f(middle);
if abs(mVal)<tolerance
    r = middle;
elseif mVal>0
    r = solveByBisection(f,lowerGuess,middle,lVal,mVal,tolerance);
else
    r = solveByBisection(f,middle,upperGuess,mVal,uVal,tolerance);
end    
end