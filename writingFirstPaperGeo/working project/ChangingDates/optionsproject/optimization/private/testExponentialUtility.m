function testExponentialUtility()
%TESTEXPONENTIALUTILITY Summary of this function goes here
%   Detailed explanation goes here
global plotsInTests

gammas = [0.5, 1.0, 2.0 ];
xValues = -1:0.01:5;
utils = zeros(length(gammas),length(xValues));
legendText = cell( length(gammas),1);
for i=1:length(gammas)
    euf = ExponentialUtilityFunction(gammas(i));
    utils(i,:) = euf.weightedUtility(xValues);
    legendText{i} = sprintf('Gamma=%f',gammas(i));
end

if plotsInTests
    figure('Name','Exponential utility functions','NumberTitle','off');
    plot( xValues, utils );
    legend( legendText, 'Location', 'SouthEast' );
end

end

