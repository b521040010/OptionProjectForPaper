function testExponentialOptimizationProblem()
    testExponentialOptimizationProblem1()
end

function testExponentialOptimizationProblem1()

    initialCash = 1;

    dayData = DayData( '20160408T150000' );

    dayData.clearInstruments();
    T = dayData.getT();
    r = dayData.getInterestRate();
    K = 1010;

    model = BlackScholesModel();
    model.S0 = 1000.0;
    model.mu = 0.08;
    model.sigma = 0.2;
    model.T = 1/12;
    callPrice = model.price(r,false,K);
    putPrice = model.price(r,true,K);

    dayData.addInstrument( Bond( T,r,1,1,Inf,Inf));
    dayData.addInstrument( Future( exp(r*T)*model.S0*0.999,exp(r*T)*model.S0*1.001, Inf, Inf));
    dayData.addInstrument( CallOption( K, callPrice*0.998, callPrice*1.002, Inf, Inf ));
    dayData.addInstrument( PutOption( K, putPrice*0.998, putPrice*1.002, Inf, Inf ));        

    assertApproxEqual( dayData.askPrice(0,DayData.futureType),exp(r*T)*model.S0*1.001,1e-6);

    currentPortfolio = OldPortfolio();
    currentPortfolio.strikes = 0;
    currentPortfolio.quantities = initialCash;
    currentPortfolio.instrumentTypes = DayData.cashType;

    problem = ExponentialOptimizationProblem(dayData, model );
    problem.currentPortfolio = currentPortfolio;

    problem = problem.computeConstants( 1 );
    problem.validateObjective(1);
    problem.validateObjective(3);
    problem.validateObjective(5);

    [~,utility1,q] = problem.optimizeFmincon();
    dq = problem.ddd.doubleQ(q);
    [utility2, ~] = problem.computeUtilityUsingMatrices(dq);
    assertApproxEqual(utility1, utility2, 1e-4);
    assert( utility1 > problem.utilityFunction.weightedUtility( 1 ));
    
    %problem.dispMatrices();
    
    [newPortfolio,utility3]=problem.optimizeExponential(initialCash);
    assert(newPortfolio.value( dayData )<initialCash);
    assert(newPortfolio.cost( dayData )>newPortfolio.value( dayData ));
    assertApproxEqual(utility1, utility3, 0.0001);
    % Note this is a pretty convincing test of the exponential optimizer
    % since we've found the same optimum using fmincon and the exponential
    % method.

    ump = problem.createUtilityMaximizationProblem();
    zeroQs = 0*q;
    zUtility1 = ump.utilityForQuantities(zeroQs);
    zUtility2 = problem.computeExpectedUtility( zeroQs');
    assertApproxEqual( zUtility1, zUtility2, 0.0001);
    
    % Now check that the solution we have computed already has the correct
    % utility and passes the constraints
    utility4 = ump.utilityForQuantities(q);
    assertApproxEqual( utility3, utility4, 0.0001);
    ump.assertConstraintsPassed(q, 0.00001);
    
    % Compute the solution using the more generic code and check
    % the answer is as we expect
    [utility5,qDash] = ump.optimize();    
    utility6 = ump.utilityForQuantities(qDash);
    assertApproxEqual(utility5, utility6, 0.001);
    assertApproxEqual(utility1, utility5, 0.001);
    ump.assertConstraintsPassed(qDash, 0.0002);
    
    % Check that our portfolio constraints are met when written
    % in an alternative manner
    dqDash = problem.ddd.doubleQ(qDash);
    newPortfolio = problem.constructNewPortfolio(dqDash);
    assert(newPortfolio.value( dayData )<initialCash); % Compare value rather than cost to cope with rounding errors
        
    % This gives yet another route to computing the optimum
end






