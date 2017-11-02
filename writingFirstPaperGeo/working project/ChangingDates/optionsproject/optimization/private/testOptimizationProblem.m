function testOptimizationProblem()
%TESTCOMPUTEEXPECTEDUTILITY Test the optimization problem class
testGeneralFunctionality();
testMosekOptimization();
%we wont use OptimizationProblem class for now
%testInvestingInCashAndStock();
%testInvestingInOptions();
end

    function testGeneralFunctionality()

        dayData = DayData( '20160408T150000' );
        dayData.clearInstruments();
        
        T = dayData.getT();
        r = dayData.getInterestRate();        
        
        model = BlackScholesModel();
        model.S0 = 1000.0;
        model.mu = 0.08;
        model.sigma = 0.1;
        model.T = 1/12;
        
        
        dayData.addInstrument( Bond( T,r,1,1,Inf,Inf));
        dayData.addInstrument( Future2( model.S0*0.9999,0,0,Inf,Inf));
        %dayData.addInstrument( Future( model.S0*0.9999,model.S0*1.000));        
        assertApproxEqual( dayData.askPrice(model.S0*0.9999,DayData.futureType),0,1e-6);

        currentPortfolio = OldPortfolio();
        currentPortfolio.strikes = 0;
        currentPortfolio.quantities = 0;
        currentPortfolio.instrumentTypes = DayData.cashType;
        
        problem = OptimizationProblem(dayData, model );
        problem.currentPortfolio = currentPortfolio;
        problem.simpleInterest = 0.01;

        assertApproxEqual( problem.computePayoff([1 0 0 0],1010), 1.01, 1e-6);
        assertApproxEqual( problem.computePayoff([0 1 0 0],1010), -1.01, 1e-6);
        %assertApproxEqual( problem.computePayoff([0 0 1 0],1010), 1010, 1e-6);
        assertApproxEqual( problem.computePayoff([0 0 1 0],1010), 1010-model.S0*0.9999, 1e-6);
        %assertApproxEqual( problem.computePayoff([0 0 0 1],1010), -1010, 1e-6);
        assertApproxEqual( problem.computePayoff([0 0 0 1],1010), model.S0*0.9999-1010, 1e-6);


        expectedUtility = problem.computeExpectedUtility( [1; 0 ]);
        payoff = (1+problem.simpleInterest)*1;
        assertApproxEqual( expectedUtility, 1-exp(-payoff), 1e-6);
        
        expectedUtility2 = problem.computeExpectedUtility( [1; 0 ], problem.getQuadRule());
        assertApproxEqual(expectedUtility,expectedUtility2,1e-4);
        

        u1 = problem.computeExpectedUtility([0 1/1000]);
        model.sigma = 0.00001;
        model.mu = 0.10;
        problem.model = model;
        u2 = problem.computeExpectedUtility([0 1/1000]);
        %assert(u1>expectedUtility);
        assert(u2>u1);

    end
    
    


    function testInvestingInCashAndStock() 
        %   Tests that only require investing in stock and cash
        %


        dayData = DayData( '20160408T150000' );
        dayData.clearInstruments();
        
        T = dayData.getT();
        r = dayData.getInterestRate();        
        
        model = BlackScholesModel();
        model.S0 = 1000.0;
        model.mu = 0.08;
        model.sigma = 0.1;
        model.T = 1/12;
        
        
        dayData.addInstrument( Bond( T,r,1,1));
        dayData.addInstrument( Future2( model.S0*0.9999,0,0)); 
        %dayData.addInstrument( Future( model.S0*0.9999,model.S0*1.000));        
        assertApproxEqual( dayData.askPrice(model.S0*0.9999,DayData.futureType),0,1e-6);
        
        currentPortfolio = OldPortfolio();
        currentPortfolio.strikes = 0;
        currentPortfolio.quantities = 1;
        currentPortfolio.instrumentTypes = DayData.cashType;
                
        problem = OptimizationProblem(dayData, model );
        problem.currentPortfolio = currentPortfolio;

        portfolio1 = problem.optimizeFmincon();
        cost = portfolio1.cost( dayData );
        assertApproxEqual(cost,1,1e-6 );
        value = portfolio1.value( dayData );
        assert( value<=cost );

        model.sigma = 0.05;
        model.mu = 0.09;
        problem.model = model;
        portfolio2 = problem.optimizeFmincon();

        assert( portfolio2.quantity(0,DayData.futureType) > ...
                portfolio1.quantity(0,DayData.futureType) );

        assert( portfolio2.quantity(0,DayData.cashType) < ...
                portfolio1.quantity(0,DayData.cashType) );

    end

    
    function testInvestingInOptions()        
    
        global plotsInTests

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
        
        dayData.addInstrument( Bond( T,r,1,1));
        dayData.addInstrument( Future( model.S0*0.999,model.S0*1.001));
        dayData.addInstrument( CallOption( K, callPrice*0.998, callPrice*1.002 ));
        dayData.addInstrument( PutOption( K, putPrice*0.998, putPrice*1.002 ));        
        assertApproxEqual( dayData.askPrice(0,DayData.futureType),1001.0,1e-6);


        currentPortfolio = OldPortfolio();
        currentPortfolio.strikes = 0;
        currentPortfolio.quantities = 1;
        currentPortfolio.instrumentTypes = DayData.cashType;
        
        problem = OptimizationProblem(dayData, model );
        problem.currentPortfolio = currentPortfolio;

        newPortfolio = problem.optimizeFmincon();
        cost = newPortfolio.cost( dayData );
        assert( cost<=1 + 1e-6 );
        value = newPortfolio.value( dayData );
        assert( value<=cost );

        if plotsInTests
            newPortfolio.plotPayoff( 'testOptimizationProblem: Optimal option portfolio', dayData.getInterestRate());
        end
        


    end

    function testMosekOptimization()

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
        %dayData.addInstrument( Future( model.S0*0.999,model.S0*1.001));
        dayData.addInstrument( Future2( model.S0*0.999,0,0,Inf,Inf));
        dayData.addInstrument( Future2( model.S0*1.001,0,0,Inf,Inf));
        dayData.addInstrument( CallOption( K, callPrice*0.999, callPrice*1.001,Inf,Inf ));
        dayData.addInstrument( PutOption( K, putPrice*0.999, putPrice*1.001, Inf, Inf ));        
        assertApproxEqual( dayData.askPrice(model.S0*0.999,DayData.futureType),0,1e-6);

        currentPortfolio = OldPortfolio();
         currentPortfolio.strikes = 0;
         currentPortfolio.quantities = 1;
         currentPortfolio.instrumentTypes = DayData.cashType;
 
         
         problem = OptimizationProblem(dayData, model );
         problem.currentPortfolio = currentPortfolio;
% 
         ump = problem.createUtilityMaximizationProblem();   
         
         ump.addConstraint(QuantityConstraint(2,-inf,0));
         ump.addConstraint(QuantityConstraint(2,0,inf));
        
%         
          expectedUtility = 1-exp(-2*(1+problem.simpleInterest));
          utilityA = ump.utilityForQuantities([1,0,0,0,0]');
          utilityB = problem.computeExpectedUtility([1,0,0,0,0]');
          assertApproxEqual(utilityA, expectedUtility, 0.001);
          assertApproxEqual(utilityA, utilityB, 0.001);
%         
%         
         [utility2,q2] = ump.optimize()
%         
%         [portfolio1,utility1, q1] = problem.optimizeFmincon([ -1;1/1000;0;0 ]);
%         cost = portfolio1.cost( dayData );
%         assert( cost<=1 + 1e-6 );
%         value = portfolio1.value( dayData );
%         assert( value<=cost );
%         
%         utility1A = ump.utilityForQuantities(q1);
%         utility1B = problem.computeExpectedUtility(q1);
%         utility2A = ump.utilityForQuantities(q2);
%         utility2B = problem.computeExpectedUtility(q2);
%         
%         ump.assertConstraintsPassed( q1, 0.005 ); %change from 0.0001
%         ump.assertConstraintsPassed( q2, 0.005 );
%                 
%         assertApproxEqual(utility1,utility2,0.002);
%         assertApproxEqual( utility2A, utility1, 0.002);
%         assertApproxEqual( utility2B, utility1, 0.002);
%         assertApproxEqual( utility1B, utility2, 0.002);
%         assertApproxEqual( utility1A, utility2, 0.002);
%         
%         portfolio2 = problem.constructNewPortfolio( buyAndSellQuantities(q2) );
%         portfolio2.disp();
%         cost = portfolio2.cost( dayData );
%         disp( cost );
%         assert( cost<=1 + 3e-3 );
%         value = portfolio2.value( dayData );
%         assert( value<=cost );                
        
    end
   
    function q= buyAndSellQuantities( netQ )
        % convert a set of net quantities into a vector of twice
        % the length containing buy and sell quantities
        q = zeros( 1, 2*length(netQ) );
        for i=1:length(netQ)
            if netQ(i)>=0
                q(2*i)=netQ(i);
            else
                q(2*i-1)=-netQ(i);
            end
        end
    end
    