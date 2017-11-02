function optionData = loadOptionData( bloomberg, fromDate, toDate, maturityDate )
% Reads all the option data from Bloomberg on a particular day for the desired
% maturity and returns
% it in a palatable format

fromDate = sprintf('%d/%d/%d',fromDate(2), fromDate(3), fromDate(1));
toDate = sprintf('%d/%d/%d',toDate(2), toDate(3), toDate(1));

fprintf('Loading data from %s-%s\n', fromDate, toDate );
data = history(bloomberg, 'SPX Index', {'PX_BID', 'PX_ASK'}, fromDate, toDate, 'Daily'  );
if (isempty(data))
    optionData = [];
    return;
end

n = size(data,1);
optionData.tickers = cell(n,1);
optionData.strikes = zeros(n,1);
optionData.bids = data(:,2);
optionData.asks = data(:,3);
optionData.isPut = zeros(n,1);
optionData.isFuture = zeros(n,1);
optionData.dates = data(:,1);

for i=1:n
    optionData.tickers{i} = 'SPX Index';
end

if mod(maturityDate(2),3)==0
    % There's a future maturing on that date
    switch maturityDate(2)
        case 3
            code = 'H';
        case 6
            code = 'M';
        case 9
            code = 'U';
        case 12
            code = 'Z';
        otherwise
            error('Invalid month %d',maturityDate(2));
    end
    number = maturityDate(1);
    if number==year(date())
        number = mod(number,10);
    else
        number = mod(number,100);
    end
    ticker = sprintf('SP%s%d Index',code,number);
    disp('Load data for ticker');
    disp(ticker);
    
    data = history(bloomberg, ticker, {'PX_BID', 'PX_ASK'}, fromDate, toDate, 'Daily'  );

    dates = data(:,1);
    bids = data(:,2);
    asks = data(:,3);
    n = size(data,1);
    tickers = cell(n,1);
    for i=1:n
        tickers{i} = ticker;
    end
    
    optionData.tickers = vertcat( optionData.tickers, tickers );
    optionData.strikes = vertcat( optionData.strikes, zeros(n,1) );
    optionData.asks = vertcat( optionData.asks, asks );
    optionData.bids = vertcat( optionData.bids, bids );
    optionData.isPut = vertcat( optionData.isPut, zeros(n,1) );
    optionData.isFuture = vertcat( optionData.isFuture, ones(n,1) );
    optionData.dates = vertcat( optionData.dates, dates );
end

expDate = sprintf('%d/%d/%d', maturityDate(2), maturityDate(3), maturityDate(1) );
price = 5*floor(median( optionData.bids )/5);
strikes = -500:5:500 + price;
chainTickers = cell([2*length(strikes), 1]);
for i=1:length( strikes )
    chainTickers{2*i-1} = sprintf('SPX US %s C%d Index', expDate, strikes(i) );
    chainTickers{2*i} = sprintf('SPX US %s P%d Index', expDate, strikes(i) );
end

[data, sec] = history( bloomberg, chainTickers, {'PX_BID', 'PX_ASK', 'PX_VOLUME', 'OPT_STRIKE_PX'}, fromDate, toDate, 'daily' );

% data and sec now contain cell arrays, to be accessed with {}
nrecords = 0;
for i=1:length(sec)   
   tickerData = data{i};   
   if ~isempty(tickerData)
       nrecords = nrecords+size(tickerData,1);
   end
end

asks = zeros( nrecords, 1 );
bids = zeros( nrecords, 1 );
strikes = zeros( nrecords, 1 );
isPut = zeros( nrecords, 1 );
isFuture = zeros( nrecords, 1 );
dates = zeros( nrecords, 1 );
tickers = cell( nrecords, 1 );

index = 1;
for i=1:length(sec)        
   tickerData = data{i};
   if ~isempty(tickerData)
       [strike, p ]=extractStrikeFromTicker( sec{i} );
       for j=1:size(tickerData,1)
           tickers{index} = sec{i};
           dates(index) = tickerData(j,1);
           bids(index) = tickerData(j,2);
           asks(index) = tickerData(j,3);
           strikes(index) = strike;
           isPut(index) = p;
           isFuture(index) = 0;
           index = index+1;
       end
   end
end    

optionData.tickers = vertcat( optionData.tickers, tickers );
optionData.strikes = vertcat( optionData.strikes, strikes );
optionData.asks = vertcat( optionData.asks, asks );
optionData.bids = vertcat( optionData.bids, bids );
optionData.isPut = vertcat( optionData.isPut, isPut );
optionData.isFuture = vertcat( optionData.isFuture, isFuture );
optionData.dates = vertcat( optionData.dates, dates );
    
assert( length(optionData.strikes)==length(optionData.asks));
assert( length(optionData.strikes)==length(optionData.bids));
assert( length(optionData.strikes)==length(optionData.isPut));
assert( length(optionData.dates)==length(optionData.dates));
    
end


