function testUtilityMaximizationProblem()
    testOptimize2();
    testUpdatePosition();
    testOptimize1();
end

function testUpdatePosition()

    prob = UtilityMaximizationProblem();
    scenarios = [90 100 110]';
    weights = [0.25 0.5 0.25]';
    prob.setQuadRule( scenarios, log(weights) );
    zcb = Bond(1,0,1,1,Inf,Inf);
    
    prob.addInstrument( zcb );
    prob.addInstrument( Future(98,99,Inf,Inf));
    %prob.setCurrentPosition( zcb );
    %%%%%%%%%%%%%%%%%
    currentPort=Portfolio();
    currentPort.add([1],{zcb})
    prob.setCurrentPosition(currentPort);
    %%%%%%%%%%%%%%%%%

    puf = PowerUtilityFunction( 1 );
    prob.setUtilityFunction( puf );    
    
    u1 = prob.optimize();

    resetData = prob.addToCurrentPosition( 1, zcb );
    u2 = prob.optimize();
    assert( u2 > u1);
    
    prob.resetPosition( resetData );
    u3 = prob.optimize();
    assertApproxEqual(u1,u3,1e-6);    

end

function testOptimize1() 
    % Basic tests that the optimize function passes

    prob = UtilityMaximizationProblem();
    scenarios = [90 100 110]';
    weights = [0.25 0.5 0.25]';
    prob.setQuadRule( scenarios, log(weights) );
    zcb = Bond(1,0,1,1,Inf,Inf);
        %%%%%%%%%%%%%%%%%
    currentPort=Portfolio();
    currentPort.add([0],{zcb})
    prob.setCurrentPosition(currentPort);
    %%%%%%%%%%%%%%%%% 
    prob.addInstrument( zcb );
    prob.addInstrument( Future(98,99,Inf,Inf));
    prob.addInstrument( CallOption(100,2,3,Inf,Inf));

    puf = PowerUtilityFunction( 3 );
    prob.setUtilityFunction( puf );    
    
    u = prob.utilityForQuantities([0 0 0]');
    assert( u==-Inf);
   
    %prob.setCurrentPosition( zcb );
    %%%%%%%%%%%%%%%%%
    currentPort=Portfolio();
    currentPort.add([1],{zcb});
    prob.setCurrentPosition(currentPort);
    %%%%%%%%%%%%%%%%% 
    lowerBound = prob.utilityForQuantities([0 0 0]')    
    %assert( lowerBound>-Inf);
    expectedLowerBound = puf.weightedUtility( 1, 0 )
    assertApproxEqual( lowerBound, expectedLowerBound, 0.0001 );
    
    qExample = [-1 1/99 0]';
    prob.assertConstraintsPassed(qExample, 0.001);    
    futureUtility = prob.utilityForQuantities(qExample);    
    expectedFutureUtility = sum(puf.weightedUtility(1/99*scenarios,log(weights)));
    assertApproxEqual( futureUtility, expectedFutureUtility, 0.0001);
    
    
    [u,q] = prob.optimize();
    assert( u >= lowerBound );
    assert( u >= futureUtility );
    prob.assertConstraintsPassed(q, 0.01);
    u2 = prob.utilityForQuantities(q);
    assertApproxEqual(u,u2,0.0001);
    
    [uDash,qDash] = prob.optimizeUnscaled();
    assertApproxEqual(u,uDash,0.0001);
    for i=1:length(q)
        assertApproxEqual(q(i),qDash(i),0.0001);
    end

    % Check that if we bound the quanties available then this
    % ceases to be an optimal solution and the new constraint is binding
    idx = 2;
    factor = 0.5;
    assert(q(idx)>0);
    prob.addConstraint( QuantityConstraint(idx,-Inf,q(idx)*factor));
    [uDash,qDash] = prob.optimizeUnscaled();
    assert(uDash<u);
    assertApproxEqual(q(idx)*factor,qDash(idx),0.0001);
    

end

function testOptimize2() 
    % Tests of optimization with power utility

    T = 1;
    r = 0.03;
    eta = 3;
    
    prob = UtilityMaximizationProblem1D();
    bsModel = BlackScholesModel();   
    bsModel.T = T;
    prob.setModel( bsModel );    
    
    utilityFunction = PowerUtilityFunction( eta );
    prob.setUtilityFunction(utilityFunction);
    
    zcb = Bond(T,r,1,1,Inf,Inf);
    %prob.setCurrentPosition( zcb );
    %%%%%%%%%%%%
    currentPort=Portfolio();
    currentPort.add([1],{zcb})
    prob.setCurrentPosition(currentPort);
    %%%%%%%%%%%%%%
    prob.addInstrument( zcb );
    
    [utility, quantities] = prob.optimize();
    assertApproxEqual( quantities(1),0,1e-6);
    x = exp(r*T)*1;
    assertApproxEqual( (x^(1-eta)-1)/(1-eta), utility, 1e-5);
    
    prob.addInstrument( Future(bsModel.S0*0.999,bsModel.S0*1.001,Inf,Inf));
    [utility2, quantities2] = prob.optimize();
    assert( utility2 > utility);
    assert( quantities2(2)>0);
    assert( quantities2(1)+1-quantities2(2)*bsModel.S0*1.001<quantities(1)+1);

end


