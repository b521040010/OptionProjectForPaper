function testBachelierModel()

S0 = 100;
meanT = 105;
sdT = 20;
T = 1;
r = 0;
gamma = 2;
principal = 0;

bm = BachelierModel();
bm.S0 = S0;
bm.meanT = meanT;
bm.sdT = sdT;

prob = UtilityMaximizationProblem1D();
prob.setModel(bm);
zcb = Bond(T,r,1,1,Inf,Inf);
%%%%%%%%%%%%%%%%%
currentPort=Portfolio();
currentPort.add([0],{zcb})
prob.setCurrentPosition(currentPort);
%%%%%%%%%%%%%%%%% 
prob.addInstrument( zcb );
prob.addInstrument( Future(exp(r*T)*S0,exp(r*T)*S0,Inf,Inf));
prob.setUtilityFunction( ExponentialUtilityFunction( gamma ));

prob.assertConstraintsPassed([0;0],1e-6);
[u,quantities] = prob.optimize();

u0 = prob.utilityForQuantities([0;0]);
u0Dash = bm.expectedExponentialUtility( principal, principal, r, T, gamma);
assertApproxEqual(u0,u0Dash,1e-6);

bondInvestment = bm.optimizeExponentialUtility( principal, r, T, gamma );
uDash = bm.expectedExponentialUtility( principal, bondInvestment, r, T, gamma);
assertApproxEqual( bondInvestment, quantities(1), 0.001);
assertApproxEqual( u, uDash, 0.001);

newPrincipal = 1;
prob = prob.createDelegate();
%%%%%%%%%%%%%%%%%
currentPort=Portfolio();
currentPort.add([newPrincipal],{zcb})
prob.setCurrentPosition(currentPort);
%%%%%%%%%%%%%%%%% 
%prob.addToCurrentPosition( newPrincipal-principal, zcb );
u0NewPrincipal = prob.utilityForQuantities([0;0]);
u0NewPrincipalDash = bm.expectedExponentialUtility( newPrincipal-principal, newPrincipal-principal, r, T, gamma);
assertApproxEqual(u0NewPrincipal,u0NewPrincipalDash,1e-6);

[u2, quantities2] = prob.optimize();

bondInvestment2 = bm.optimizeExponentialUtility( newPrincipal, r, T, gamma );
uDash2 = bm.expectedExponentialUtility( newPrincipal, bondInvestment2, r, T, gamma);
assertApproxEqual( bondInvestment2, newPrincipal-quantities(2)*exp(r*T)*S0, 0.01);
assertApproxEqual( u2, uDash2, 0.001);

assert( u2>u); % If you have more money to start with, you're happier at the end

end

