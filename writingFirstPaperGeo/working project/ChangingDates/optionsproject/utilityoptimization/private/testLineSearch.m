function testLineSearch()
%TESTLINESEARCH Summary of this function goes here
%   Detailed explanation goes here

f = @(x) x^2 - 2;

root2 = onlyLineSearch(f,1.35,1e-6);
assertApproxEqual(root2,sqrt(2),1e-6);



end

