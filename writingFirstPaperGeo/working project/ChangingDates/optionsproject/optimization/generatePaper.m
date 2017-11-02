function generatePaper()
%GENERATEPAPER Creates the plots required for the paper on
%  indifference pricing
initutilityoptimization();
delete(findall(0,'Type','figure'));
mkdir('plots');

lastPlot = plotDeltaHedgeIndifferencePrice();
saveFigure(lastPlot,'plots/indifferencePriceDeltaHedge.pdf', [20 15]);

lastPlot = plotCallIndifferencePrice( createParams(), chooseStrike());
saveFigure(lastPlot,'plots/indifferencePrice.pdf', [20 15]);

function [newParams, title] = riskSensitivity( params, rho )
    newParams = params;
    newParams.riskAversion  = rho;
    title = sprintf('Risk aversion %0.1f', rho);
end
lastPlot = plotIndifferencePriceSensitivities( createParams(), chooseStrike(), @riskSensitivity, [100 4 2 1 0.5 0.3 0.1] );
saveFigure(lastPlot,'plots/indifferencePriceByRiskAversion.pdf', [20 15]);


function [newParams, title] = modelSensitivity( params, sigmaMultiplier )
    newParams = params;
    newParams.riskAversion  = sigmaMultiplier;
    title = sprintf('Volatility multiplier %0.1f', sigmaMultiplier);
end
lastPlot = plotIndifferencePriceSensitivities( createParams(), chooseStrike(), @modelSensitivity, [0.25 0.5 1 2.0 4.0] );
saveFigure(lastPlot,'plots/indifferencePriceByVolatility.pdf', [20 15]);


function [newParams, title] = positionSensitivity( params, callQuantity )
    newParams = params;
    newParams.digitalCallQuantity  = callQuantity;
    title = sprintf('%0.1f digital calls held', callQuantity);
end
lastPlot = plotIndifferencePriceSensitivities( createParams(), chooseStrike(), @positionSensitivity, [-0.2 -0.1 0 0.1 0.2] );
saveFigure(lastPlot,'plots/indifferencePriceByPosition.pdf', [20 15]);

params = createParams();
[plotA,plotB] = plotOptimalPortfolio(params);
disp('A mildy risk averse investor takes a short strangle');
saveFigure(plotA,'plots/blackScholesOptimalPortfolioA.pdf', [20 15]);
saveFigure(plotB,'plots/blackScholesOptimalPortfolioB.pdf', [20 15]);

params = createParams();
params.model = StudentTModel();
plotA = plotOptimalPortfolio(params);
disp('With fat tails we take a more complex position');
saveFigure(plotA,'plots/studentOptimalPortfolio.pdf', [20 15]);

params = createParams();
params.model = StudentTModel();
params.riskAversion = 5.0;
plotA = plotOptimalPortfolio(params);
disp('As risk aversion increases, we hedge the tails more carefully');
saveFigure(plotA,'plots/studentOptimalPortfolioRiskAverse.pdf', [20 15]);

params = createParams();
params.model = StudentTModel();
plotA = plotOptimalPortfolioNoTails(params, 1.0);
saveFigure(plotA,'plots/studentOptimalPortfolioHistoric.pdf', [10 10]);
disp('If we believe historic volatility');
plotA = plotOptimalPortfolioNoTails(params, 3.0);
saveFigure(plotA,'plots/studentOptimalPortfolioHighVol.pdf', [10 10]);
disp('If we believe volatility will grow');
plotA = plotOptimalPortfolioNoTails(params, 1./3.0);
saveFigure(plotA,'plots/studentOptimalPortfolioLowVol.pdf', [10 10]);
disp('If we believe volatility will shrink');

end

function [proportions, filter] = choosePlotProportions()
% Standard plot points and standard positions for markers
    proportions = [-0.5, -0.4, -0.3, -0.2, -0.1, -0.01, 0,  0.01, 0.1, 0.2, 0.3, 0.4, 0.5];
    filter =      [0      0     1     0     1    0      0   0     0    1    0    1    0    ];
    filter = filter==1;
end

function strike = chooseStrike() 
% Standard strike
    strike = 1990;
end

function params = createParams()
% Create default parameters to pass to other functions
    params.riskAversion = 1.0;
    params.date = '2014-09-09';
    params.K = chooseStrike();
    params.kRange = 200;
    params.model = BlackScholesModel();
    params.ewmaLambda = 0.95;
    params.sigmaMultiplier = 1.0;
    params.digitalCallQuantity = 0.0;
end

function [plotA, plotB] = plotOptimalPortfolio( params )
% Plot the optimal portfolio
    ump = createProblem(params);
    addQuantityConstraints(ump);
    [utility, quantities] = ump.optimize();
    plotA = ump.plotPortfolio( sprintf('Optimal portfolio, BS model. Risk Aversion %0.1f, Utility=%0.2f', params.riskAversion, utility),quantities);
    plotB = plotPortfolio( ump.getInstruments(), quantities);
end

function [plotA,plotB] = plotOptimalPortfolioNoTails( params, sigmaMultiplier)
% Plot the optimal portfolio ignoring far out of the money options
    params.sigmaMultiplier = sigmaMultiplier;
    [ump, ~] = createProblem(params);
    removeOutOfTheMoney(ump,params.K,params.kRange);
    addQuantityConstraints(ump);
    [utility, quantities] = ump.optimize();
    plotA = ump.plotPortfolio( sprintf('Sigma multiplier %0.1f, Utility=%0.2f', sigmaMultiplier, utility),quantities, [-7 7]);
    plotB = plotPortfolio( ump.getInstruments(), quantities);
end

function plotNum = plotCallIndifferencePrice( params, strike )

    [ump, dayData] = createProblem(params);
    removeOutOfTheMoney(ump,params.K,params.kRange);

    call = dayData.findInstrument( strike, DayData.callType );
    put = dayData.findInstrument( strike, DayData.putType );
    ump.removeInstrument( call );
    ump.removeInstrument( put );

    addQuantityConstraints(ump);

    proportions = choosePlotProportions();
    indifferencePrices = zeros( 1, length( proportions ));
    zcb = dayData.findInstrument(0, DayData.cashType );
    for i=1:length( proportions)
        price = call.getAsk();
        proportion = proportions(i);
        indifferencePrices(i) = ump.indifferencePriceExponentialUtility(zcb, proportion/price, call, proportion);
    end
    

    plotNum = figure();
    subplot(2,1,2);
    % We use seller's prices for this report, hence the minuses
    plot( -proportions, indifferencePrices*price./proportions, '-b' );
    hold on;
    plot( -proportions, call.getBid()*ones( size( proportions )), '--r');
    plot( -proportions, call.getAsk()*ones( size( proportions )), '--k');
    xlabel('Proportion of $1 ask quantity');
    ylabel('Price per unit');
    legend('Indifference price', 'Bid', 'Ask', 'Location','northwest');
    title( sprintf('Indifference price of call option with strike %d', strike));
    hold off;
    subplot(2,1,1);
    plot( -proportions, -indifferencePrices, '-b');
    hold on;
    plot( -proportions, -proportions*call.getBid()/price, '--r');
    plot( -proportions, -proportions*call.getAsk()/price, '--k');
    xlabel('Proportion of $1 ask quantity');
    ylabel('Price of contract');
    legend('Indifference price', 'Bid', 'Ask', 'Location','northwest');
    title( sprintf('Indifference price of call option with strike %d', strike));
    hold off;
    
end

function plotNum = plotIndifferencePriceSensitivities( params, strike, sensitivityFunction, sensitivities )

    sensitivities = sort( sensitivities );
    [ump, dayData] = createProblem(params);

    call = dayData.findInstrument( strike, DayData.callType );
    put = dayData.findInstrument( strike, DayData.putType );

    [proportions, filter] = choosePlotProportions();
    
    plotNum = figure();
    l{1} = 'Bid';
    l{2} = 'Ask';
    l{3} = 'Subhedge';
    l{4} = 'Superhedge';
    plot( proportions, call.getBid()*ones( size( proportions )), '--r');
    hold all;
    plot( proportions, call.getAsk()*ones( size( proportions )), '--k');    
    
    removeOutOfTheMoney(ump,params.K,params.kRange);
    ump.removeInstrument( call );
    ump.removeInstrument( put );
    [sup,sub] = ump.superHedgePrice( call );
    plot( proportions, sub*ones( size( proportions )), ':r');    
    plot( proportions, sup*ones( size( proportions )), ':k');    
    
    cmap = colormap('cool');
    
    % We have to do two plots, one for markers one for lines
    count = 5;
    for pass=1:2
        if (pass==1)
            filteredProportions = proportions(filter);
        else
            filteredProportions = proportions;
        end
        for idx = 1:length(sensitivities)
            s = sensitivities(idx);

            [newParams,label] = sensitivityFunction( params, s );
            ump = createProblem(newParams);
            removeOutOfTheMoney(ump,params.K,params.kRange);

            ump.removeInstrument( call );
            ump.removeInstrument( put );

            addQuantityConstraints(ump);

            indifferencePrices = zeros( 1, length( filteredProportions ));
            zcb = dayData.findInstrument(0, DayData.cashType );
            for i=1:length( filteredProportions)
                price = call.getAsk();
                proportion = filteredProportions(i);
                indifferencePrices(i) = ump.indifferencePriceExponentialUtility(zcb, proportion/price, call, proportion);
            end

            lambda = (idx-1)/(length(sensitivities)-1);
            cIdx = round((size(cmap,1)-1)*lambda+1);
            color = cmap(cIdx,:);
            if lambda<0.49
                marker = 'v';
            elseif lambda>0.51
                marker = '^';
            else
                marker = 'o';
            end
            
            markerSize = (abs(lambda-0.5)^2)*20 + 5;            
            
            linespec = '-'; 
            if (pass==1)
                linespec = marker;
                l{count} = label;
                count = count+1;
            end
            % We use sellers prices hence the minus signs below
            plot( -filteredProportions, indifferencePrices*price./filteredProportions, linespec, 'Color', color, 'MarkerSize', markerSize);

        end
    end

    xlabel('Proportion of $1 ask quantity');
    ylabel('Price per unit');
    legend(l,'Location','northwest');
    title( sprintf('Indifference price of call option with strike %d', strike));
    hold off;
end



function [ump, dayData, model] = createProblem( params )
% Create a basic problem instance associated with the given date
% and risk aversion. The model is calibrated using the Black Scholes Model
% with EWMA

    utilityFunction = ExponentialUtilityFunction( params.riskAversion );

    dayData = DayData( params.date );
    arbitrage = ArbitrageFinder.findArbitrageForDate( params.date, false, true );
    assert(~arbitrage);

    model = dayData.calibrateHistoric( params.model, params.ewmaLambda);
    model.sigma = params.sigmaMultiplier * model.sigma;
    ump = UtilityMaximizationProblem1D();
    ump.setModel( model );
    
    % Create a utility maximization problem corresponding
    % to this problem
    ump.setUtilityFunction( utilityFunction );
    
    digitalCall = DigitalCallOption( params.K, 0, Inf);
    p = Portfolio();
    q(1) = params.digitalCallQuantity;
    i{1} = digitalCall;
    p.add( q, i );
    ump.setCurrentPosition( p );
    
    for i=1:length(dayData.instruments)
        ump.addInstrument( dayData.instruments{i} );
    end
end


function removeOutOfTheMoney(ump, K, kRange) 
    % Get rid of instruments which are too far out of the money
    function accept = filter( instrument )
        if (isa(instrument,'CallOption') || isa(instrument,'PutOption'))
            accept = abs( instrument.K - K) < kRange;
        else 
            accept = 1;
        end
    end
    ump.filterInstruments( @filter );
end

function addQuantityConstraints(ump) 
% Add a constraint that a maximum of $1 can be bought or sold of
% each instrument
    instruments = ump.getInstruments();
    nInstruments = length(instruments);
    for i=1:nInstruments
        instrument = instruments{i};
        min = -1;
        max = 1;
        if (isfinite(instrument.getBid()))
            min = -1/instrument.getBid();
        end
        if (isfinite(instrument.getAsk()))
            max = 1/instrument.getAsk();
        end
        ump.addConstraint( QuantityConstraint(i, min, max));
    end
end

function plotNum = plotDeltaHedgeIndifferencePrice()

    params = createParams();
    [ump, dayData] = createProblem(params);
    removeOutOfTheMoney(ump,params.K,params.kRange);
    strike = chooseStrike();
    future = dayData.findInstrument( 0, DayData.futureType );
    call = dayData.findInstrument( strike, DayData.callType );
    put = dayData.findInstrument( strike, DayData.putType );
    ump.removeInstrument( call );
    ump.removeInstrument( put );
    
    [proportions, filter] = choosePlotProportions();
    plotNum = figure();
    l{1} = 'Bid';
    l{2} = 'Ask';
    l{3} = 'Subhedge';
    l{4} = 'Superhedge';
    l{5} = 'Indifference Price';
    plot( proportions, call.getBid()*ones( size( proportions )), '--r');
    hold all;
    plot( proportions, call.getAsk()*ones( size( proportions )), '--k');    
    
    [sup,sub] = ump.superHedgePrice( call );
    plot( proportions, sub*ones( size( proportions )), ':r');    
    plot( proportions, sup*ones( size( proportions )), ':k');    
    
    riskAversions = [0.5 1 2];
    indifferencePrices = zeros( length(riskAversions), length( proportions ));
    for outerIdx = 1:length( riskAversions)
        riskAversion = riskAversions(outerIdx);
        for idx=1:length(proportions)
            p = proportions(idx);
            r = dayData.getInterestRate();
            bsm =dayData.blackScholesModel();
            bsm.sigma = bsm.impliedVolatility( r, strike, false, 0.5*(call.getBid()+call.getAsk()));
            dhp = DeltaHedgeProblem();
            dhp.bsm = bsm;           
            dhp.proportionTransactionCosts = (future.getAsk()-future.getBid())/future.getAsk();
            dhp.nPaths = 10000;
            dhp.K = strike;
            dhp.isPut = false;
            dhp.r = r;
            dhp.utilityFunction = ExponentialUtilityFunction( riskAversion);
            quantity = p / call.getAsk();
            priceGuess = bsm.price( dhp.r, dhp.isPut, dhp.K );
            for nsteps=1:10
                rng('default');
                S = dhp.simulatePricePaths( nsteps );
                priceForQuantity = dhp.sellersIndifferencePrice( quantity, priceGuess, S );
                if nsteps==1
                    bestPriceForQuantity = priceForQuantity;
                else                
                    %if (priceForQuantity<bestPriceForQuantity)
                    %    fprintf('Min at %d',nsteps);
                    %end
                    bestPriceForQuantity = min(priceForQuantity,bestPriceForQuantity);
                end
            end
            indifferencePrices(outerIdx, idx) = bestPriceForQuantity/quantity;
        end
    end

    % Plot the results prettily
    cmap = colormap('cool');
    count = 5;
    for pass=1:2
        for idx = 1:length( riskAversions)
            riskAversion = riskAversions(idx);
            lambda = (idx-1)/(length(riskAversions)-1);
            cIdx = round((size(cmap,1)-1)*lambda+1);
            color = cmap(cIdx,:);
            if lambda<0.49
                marker = 'v';
            elseif lambda>0.51
                marker = '^';
            else
                marker = 'o';
            end

            markerSize = (abs(lambda-0.5)^2)*20 + 5;            

            linespec = '-'; 
            if pass==1
                linespec = marker;
                l{count} = sprintf('Risk aversion %0.1f', riskAversion);
                count = count+1;
                f = filter;
            else
                f = proportions==proportions;
            end        
            plot( proportions(f), indifferencePrices(idx,f), linespec,'Color', color, 'MarkerSize', markerSize );
        end
    end

    legend(l,'Location', 'northwest');
    title( 'Indifference price using discrete delta hedge strategy');
    hold off;
    
end

   