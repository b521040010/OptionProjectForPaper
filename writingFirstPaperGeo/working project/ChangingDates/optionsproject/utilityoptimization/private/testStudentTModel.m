function testStudentTModel()
%TESTSTUDENTTMODEL Summary of this function goes here
%   Detailed explanation goes here
S0 = 100.0;
mu = log(120.0);
T = 1;
sigma = 0.2;
nu = 2.2;
tDistributed = trnd(nu,10000,1);
samplePrices = exp(tDistributed*sigma + mu);
sampleReturns = (samplePrices-S0)/S0;

model = StudentTModel();
model = model.fit( S0, T, sampleReturns, ones(size(samplePrices)));
assertApproxEqual( model.mu, mu, 0.1 );
assertApproxEqual( model.sigma, sigma, 0.02 );
assertApproxEqual( model.nu, nu, 0.5 );
assertApproxEqual( model.T, T, 0.001 );

totalProb = integral(@(x) model.pdf(x),0, Inf);
assertApproxEqual( totalProb, 1.0, 0.0001);
totalProb2 = integral(@(x) exp(model.logPdf(x)),0, Inf);
assertApproxEqual( totalProb2, 1.0, 0.0001);

wayPoints = model.getWayPoints();
quadRule = QuadRule.adapted(@(x) model.pdf(x), wayPoints);
totalProb3 = quadRule.integrate(@(x) model.pdf(x) );
assertApproxEqual( totalProb3, 1.0, 0.02);

end

