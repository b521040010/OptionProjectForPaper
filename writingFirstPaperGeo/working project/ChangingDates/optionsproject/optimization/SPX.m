classdef SPX
    %SPX Functions to work with SPX index options
    
    properties
    end
    
    methods (Static) 
        function dateVector = findExpiryDate( year, month )
            if (year<2015)
                % S & P Options expire on the "Saturday following the 3rd Friday of the
                % month". Which is not quite the same thing as the 3rd Saturday of the
                % month!
                error('SPX changed their methodology in 2015 - options now expire on Fridays');
            end
            
            cal = calendar( year, month );
            count = 0;
            i = 1;
            while (count<3)
                if cal(i,6)>0
                    count = count+1;
                    day = cal(i,6);
                end
                i = i+1;    
            end
            dateVector = [ year month day 0 0 0 ];
        end
        
        function [ ret ] = findNextExpiryDate( date, quarterly )
        %FINDNEXTEXPIRYDATE Find the next S&P option expiry date.

            if (nargin<2)
                quarterly = false;
            end
            dateVec = date;

            ret = SPX.findExpiryDate( dateVec(1), dateVec(2) );
            if (ret(3)<dateVec(3) || (quarterly && mod(dateVec(2),3)~=0))
                if quarterly
                    nextMonth = (fix(ret(2)/3) + 1)*3;
                else 
                    nextMonth = ret(2)+1;
                end
                nextYear = ret(1);
                if nextMonth>12 
                    nextMonth = nextMonth-12;
                    nextYear = nextYear +1;
                end
                ret = SPX.findExpiryDate( nextYear, nextMonth );     
            end

        end        
    end
    
end

