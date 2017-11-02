function writeOptionData( data, directory )
%WRITEOPTIONDATA Write option data that has been read from Bloomberg to a
%file

dates = unique(data.dates);
for i=1:length(dates)
    date = dates(i);
    dateStr = datestr( date, 'yyyy-mm-dd');
    fprintf('Recording data for %s\n',dateStr);
    fileName = sprintf('%s/%s.csv',directory,dateStr);
    file = fopen( fileName, 'w' );
    fprintf(file,'Ticker,Strike,Bid,Ask,isPut,isFuture\n');
    for j=1:size(data.strikes,1)
        if data.dates(j)==date
            fprintf(file,'%s,%d,%d,%d,%i,%i\n',data.tickers{j},data.strikes(j),data.bids(j),data.asks(j),data.isPut(j),data.isFuture(j));
        end
    end
    fclose( file );
end

end

