function testQuadRule()

gaussRule = QuadRule.gaussLegendre(7);
value = gaussRule.integrate(@exp);
assertApproxEqual( value, exp(1)-exp(-1),0.0001);

gaussRule1 = QuadRule.gaussLegendre(7,2,5);
value = gaussRule1.integrate(@exp);
assertApproxEqual( value, exp(5)-exp(2),0.0001);

gaussRule2 = QuadRule.gaussLegendre(7,-2,2);
rule = QuadRule.combine( [gaussRule1 gaussRule2] );
value = rule.integrate(@exp);
assertApproxEqual( value, exp(5)-exp(-2),0.0000001);

model = BlackScholesModel();
model.sigma = 0.001;
quadRule = QuadRule.adapted( @(x) x.*model.pdf(x), model.S0, model.getWayPoints() );
quadMean = quadRule.integrate( @(x) x.*model.pdf(x) );
assertApproxEqual( model.mean(), quadMean, 1e-6 );

end

