function transferHistoricDataToFile( bloomberg, folder, startDate, endDate )
% Run the command
%  transferHistoricDataToFile( bloomberg, 'SPX Index', [2010 1 1], [2013 12 1]);
%  to download lots of historic option data to a file. 


mkdir( folder );
yStart = startDate(1);
mStart = startDate(2);
dStart = startDate(3);
yEnd = endDate(1);
mEnd = endDate(2);
dEnd = endDate(3);
while (yStart<yEnd || (yStart==yEnd && mStart<mEnd) || (yStart==yEnd && mStart==mEnd && dStart<dEnd))
    date = [yStart, mStart, dStart, 0, 0, 0];
    expiryDate = findNextExpiryDate( date, true );
    optionData = loadOptionData( bloomberg, date, expiryDate, expiryDate);
    if (~isempty( optionData))
        writeOptionData( optionData, folder );
    end
    yStart = expiryDate(1);
    mStart = expiryDate(2);
    dStart = expiryDate(3)+1;
end

end

