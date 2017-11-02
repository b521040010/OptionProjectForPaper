classdef DoubledDayData
    %DOUBLEDDAYDATA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nInstruments;
        instrumentPrice;
        instrumentStrike;        
        instrumentType;   
        isShort;
    end
    
    methods
        function ddd=DoubledDayData( dd )
            n = length(dd.instruments);
            ddd.nInstruments = n*2;
            ddd.instrumentPrice = zeros(1,n*2);
            ddd.instrumentStrike = zeros(1,n*2);
            ddd.instrumentType = zeros(1,n*2);
            ddd.isShort=zeros(1,n*2);
            for i=1:n
                ins = dd.instruments{i};
                K = 0;
                type = DayData.cashType;
                if isa(ins,'Bond')
                    type = DayData.cashType;
                elseif isa(ins,'Future2')
                    type = DayData.futureType;
                    K = ins.getStrike();
                elseif isa(ins,'PutOption')
                    type = DayData.putType;
                    K = ins.getStrike();
                elseif isa(ins,'CallOption')
                    type = DayData.callType;
                    K = ins.getStrike();
                end                    
                ddd.instrumentPrice(2*i-1) = ins.getAsk();
                ddd.instrumentPrice(2*i) = -ins.getBid();
                ddd.isShort(2*i-1)=0;
                ddd.isShort(2*i)=1;
                ddd.instrumentStrike(2*i-1) = K;
                ddd.instrumentStrike(2*i) = K;
                ddd.instrumentType(2*i) = type;
                ddd.instrumentType(2*i-1) = type;
            end
        end
        
        function dq = doubleQ( ddd, netQ )
            % Given a set of net quantities double the array size
            % to get a positive vector of buys and sells
            dq = zeros(1,2*length(netQ));
            assert( length(netQ)==ddd.nInstruments/2);
            for i=1:length(netQ)
                dq(2*i-1) = (netQ(i)>0) .* netQ(i);
                dq(2*i) = (netQ(i)<0) .* -netQ(i);
            end
        end

        function netQ = netQ( ddd, dq )
            % inverse of doubleQ
            netQ = zeros(1,ddd.nInstruments/2);
            for i=1:ddd.nInstruments
                if (mod(i,2)==1)
                    netQ((i+1)/2) = netQ((i+1)/2)+dq(i);
                else
                    netQ(i/2) = netQ(i/2)-dq(i);
                end
            end
        end
        
    end


    
end

