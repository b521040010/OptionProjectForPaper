function testOnlyOptimizationNStepFromDate(fromDate,nStep)
%testOnlyOptimizationNStepFromDate('2016-04-08 15:00:00',10)
for i=1:nStep
disp(fromDate)
fromDateNum=datenum(fromDate);
fromDateNum10SecBefore=datestr(fromDateNum-10/(24*60*60),'yyyy-mm-dd HH:MM:SS');
fromDate=fromDateNum10SecBefore;
date=strcat(fromDateNum10SecBefore(1:4),fromDateNum10SecBefore(6:7),fromDateNum10SecBefore(9:10),'T',fromDateNum10SecBefore(12:13),fromDateNum10SecBefore(15:16),fromDateNum10SecBefore(18:19));
testOnlyOptimizationForDate(date);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
end
end