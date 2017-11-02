function [bid ask bidSize askSize] = loadBidAskBidsizeAsksizeGivienDateAndTime(bloomberg,ticker,dateTime)
%we want the prices 1 second before the given time unless we will not be able to
%trade


 if mod(datenum(dateTime)-1/(24*60*60),1)==0
     oneScecondBeforeDateTime=strcat(datestr(datenum(dateTime)-1/(24*60*60)),' 00:00:00');
 else
     oneScecondBeforeDateTime=datestr(datenum(dateTime)-1/(24*60*60));
 end

% %we look back timeLookBack mins before the given date
 timeLookBack=0.2;
 startdate=datestr(datenum(oneScecondBeforeDateTime)-timeLookBack*1/(24*60));
 data = timeseries(bloomberg,ticker,{startdate,oneScecondBeforeDateTime},[],{'BID','ASK'});
 
 temp=1;
 i=1;
 if cell2mat(data(size(data,1),1)) == 'BID'
     bid=cell2mat(data(size(data,1),3));
     bidSize=cell2mat(data(size(data,1),4));
     while temp==1
         if cell2mat(data(size(data,1)-i,1)) == 'ASK'
             ask=cell2mat(data(size(data,1)-i,3));
             askSize=cell2mat(data(size(data,1)-i,4));
             temp=0;
         else
             i=i+1;
         end
     end
     
 else if cell2mat(data(size(data,1),1)) == 'ASK'
     ask=cell2mat(data(size(data,1),3));
     askSize=cell2mat(data(size(data,1),4));
     while temp==1
         if cell2mat(data(size(data,1)-i,1)) == 'BID'
             bid=cell2mat(data(size(data,1)-i,3));
             bidSize=cell2mat(data(size(data,1)-i,4));
             temp=0;
         else
             i=i+1;
         end
     end
     
 end
     

end