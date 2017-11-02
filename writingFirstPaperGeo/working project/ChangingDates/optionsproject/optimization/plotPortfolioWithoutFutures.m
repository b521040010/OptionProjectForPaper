function [ fig ] = plotPortfolioWithoutFutures( date,instruments, quantities )
%PLOTPORTFOLIO Plot the quantities of call options bought and sold in 
%   a readable manner
j=1;
nInstruments = length( instruments );
nOptions = 0;
strikes = zeros(1,nInstruments);
for i=1:nInstruments
    instrument = instruments{i};
    if (isa(instrument,'CallOption') || isa(instrument,'PutOption'))
        K = instrument.K;
        nOptions = nOptions+1;
        strikes(nOptions)=K;
    end
end
strikes = sort( unique( strikes(:,1:nOptions) ));
call = zeros(1,length(strikes));
put = zeros(1,length(strikes));
askSize = zeros(1,length(strikes));
bidSize = zeros(1,length(strikes));
askSizePut = zeros(1,length(strikes));
bidSizePut = zeros(1,length(strikes));
askSizeFuture = 0;
bidSizeFuture = 0;
future = 0;
bond = 0;
for i=1:nInstruments
    instrument = instruments{i};
    if (isa(instrument,'CallOption'))
        K = instrument.K;
        index = find(strikes==K);
        call(index) = call(index) + quantities(i);
        askSize(index)=askSize(index) + instrument.askSize;
        bidSize(index)=bidSize(index) - instrument.bidSize;
    elseif (isa(instrument,'PutOption'))
        K = instrument.K;
        index = find(strikes==K);
        put(index) = put(index) + quantities(i);
        askSizePut(index)=askSizePut(index) + instrument.askSize;
        bidSizePut(index)=bidSizePut(index) - instrument.bidSize;
    elseif (isa(instrument,'Future2'))
        %future = future + quantities(i);
        future(j) = quantities(i);
        strikeFuture(j)=instrument.K;
        askSizeFuture=askSizeFuture+instrument.askSize;
        bidSizeFuture=bidSizeFuture- instrument.bidSize;
        j=j+1;
    elseif (isa(instrument,'Bond'))
        bond = bond + quantities(i);
    end
end
fig = figure();
subplot(3,1,1);
% plot(strikes,bidSize,'r');
% hold on;
% plot(strikes,askSize,'r');
% hold on;

bar(strikes, call  );
%ylim([-120 30])
xlabel('Strikes');
%ylabel('Quantity (contracts)');

title(strcat('Calls',date));
subplot(3,1,2);
% plot(strikes,bidSizePut,'r');
% hold on;
% plot(strikes,askSizePut,'r');
% hold on;
bar( strikes, put  );
xlabel('Strikes');
ylabel('Quantity (contracts)');
title('Puts');
% subplot(3,1,3);
% 
% bar( [future ;0 0], 'grouped');
% ylabel('Quantity');
% strike1=num2str(strikeFuture(1));
% strike2=num2str(strikeFuture(2));
% bidSizeFutureVec=[bidSizeFuture bidSizeFuture];
% askSizeFutureVec=[askSizeFuture askSizeFuture];
% plot(strikeFuture,bidSizeFutureVec,'r');
% hold on;
% plot(strikeFuture,askSizeFutureVec,'r');
% hold on;
% strike1=strcat('Future with Strike= ',strike1);
% strike2=strcat('Future with Strike= ',strike2);
% bar(strikeFuture, future,0.25  );
% xlabel('Strikes');
% %ylabel('Quantity (contracts)');
% title('Futures');
subplot(3,1,3); 
bar( [bond;0], 'grouped');
ylabel('Quantity');
title('Cash');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot each quantity searately
% fig = figure();
% plot(strikes,bidSize,'r');
% hold on;
% plot(strikes,askSize,'r');
% hold on;
% bar(strikes, call  );
% xlabel('Strikes');
% ylabel('Quantity');
% title('Calls');
% 
% fig = figure();
% plot(strikes,bidSizePut,'r');
% hold on;
% plot(strikes,askSizePut,'r');
% hold on;
% bar(strikes, put  );
% xlabel('Strikes');
% ylabel('Quantity');
% title('Puts');
% 
% 
% 
% fig = figure();
% strikeFuture
% bidSizeFutureVec=[bidSizeFuture bidSizeFuture]
% askSizeFutureVec=[askSizeFuture askSizeFuture]
% plot(strikeFuture,bidSizeFutureVec,'r');
% hold on;
% plot(strikeFuture,askSizeFutureVec,'r');
% hold on;
% bar(strikeFuture, future,0.25  );
% xlabel('Strikes');
% ylabel('Quantity');
% title('Futures');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end