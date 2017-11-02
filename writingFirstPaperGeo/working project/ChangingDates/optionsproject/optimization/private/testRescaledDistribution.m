function testRescaledDistribution()
%TESTRESCALEDDISTRIBUTION Summary of this function goes here
%   Detailed explanation goes here
d = rescaledDistribution();
d.mu = 5;
d.sigma = 4;
d.delegate = tDistribution();
d.delegate.nu = 3;

assertApproxEqual( d.integrate(1,0,-Inf,Inf), 5.0, 0.00001 );
assertApproxEqual( d.integrate(0,1,-Inf,Inf), 1.0, 0.00001 );

end

