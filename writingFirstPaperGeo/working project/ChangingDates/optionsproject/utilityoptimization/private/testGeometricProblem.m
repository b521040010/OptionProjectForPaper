function testGeometricProblem()
%TESTSEPARABLEPROBLEM Summary of this function goes here
%   Detailed explanation goes here

% Minimize 4 e^3x subject to 2<x<5
prob = GeometricProblem(1);
prob.setObjective( 4, 3);
prob.addLinearLowerBound(0.5,1,'Lower bound test');
prob.addLinearUpperBound(1,5,'Upper bound test');
assertApproxEqual( exp(prob.computeObjective(2)), 4*exp(3*2),0.0001);
[objective, x, ~] = prob.optimize();
assertApproxEqual(x,2,0.0001);
assertApproxEqual(exp(objective),4*exp(3*2), 0.01);

prob.assertConstraintsPassed(x, 0.001);

% Minimize e^x_1 + e^{-x_2}
prob = GeometricProblem(2);
prob.setObjective( [1 1]', [1 0; 0 -1]);
prob.addLinearLowerBound([1 0],-1, 'Lower bound test');
prob.addLinearUpperBound([0 1],1, 'Upper bound test');
assertApproxEqual( exp(prob.computeObjective([0 0]')), 2.0,0.0001);
[objective, x, ~] = prob.optimize();
assertApproxEqual(2.0*exp(-1),exp(objective),0.0001);
assertApproxEqual(x(1),-1.0, 0.01);
assertApproxEqual(x(2),1.0, 0.01);

prob.assertConstraintsPassed(x, 0.001);


end

