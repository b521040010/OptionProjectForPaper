function testLegendrePolynomial()
%TESTLEGENDREPOLYNOMIAL Summary of this function goes here
%   Detailed explanation goes here

lp3 = legendrePolynomial(3);
assertApproxEqual( lp3(1), 0, 0.0001);
assertApproxEqual( lp3(2), -3/2, 0.0001);
assertApproxEqual( lp3(3), 0, 0.0001);
assertApproxEqual( lp3(4), 5/2, 0.0001);

end

