function testBlackScholesModel()

m = BlackScholesModel();
m.S0 = 100;
assertApproxEqual( integral( @m.pdf,0,Inf), 1.0, 0.0001);
assertApproxEqual( integral( @(x) m.pdf(x),0,Inf), 1.0, 0.0001);
numericMean = integral( @(x) x.*m.pdf(x),0,Inf);
assertApproxEqual( numericMean, m.S0 * exp(m.mu * m.T), 0.0001);
assertApproxEqual( numericMean, m.mean(), 0.0001);
numericMean2 = integral( @(x) x.*exp( m.logPdf(x)),0,Inf);
assertApproxEqual( numericMean, numericMean2, 0.0001);
eXSquared = integral( @(x) x.^2.*m.pdf(x),0,Inf);
numericVariance = eXSquared - numericMean^2;
assertApproxEqual( sqrt(numericVariance), m.sd(), 0.0001);

m.T = 0.25;
m.sigma = 0.1;
m.S0 = 90;
r = 0.05;
isPut = false;
K = 100;

[bsPrice, delta] = m.price(r,isPut,K);
assert( abs(bsPrice-0.058)<0.001 );
n = m;
h = 0.001;
n.S0 = n.S0 + h;
bsPrice2 = n.price(r,isPut,K);
deltaEstimate = (bsPrice2 - bsPrice)/h;
assertApproxEqual( deltaEstimate, delta, 0.0001);


logNormalM = 0.08*m.T;
logNormalS = 0.2*sqrt(m.T);
otherModel = m.fitLogNormalReturn( logNormalM, logNormalS );
[actualM, actualS] = otherModel.logNormalParameters();
assertApproxEqual( actualS, logNormalS, 0.00001);
assertApproxEqual( actualM-log(m.S0), logNormalM, 0.00001);

rng('default');
[paths, ~] = m.simulatePricePaths(100000,5);
finalPrice = paths(:,end);
assertApproxEqual(mean( finalPrice ), m.mean(),0.1);

end

