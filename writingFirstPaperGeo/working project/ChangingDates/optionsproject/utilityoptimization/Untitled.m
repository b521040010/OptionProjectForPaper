ump = UtilityMaximizationProblem();
    zcb = Bond(1,0.1,1,1,Inf,Inf);
    currentPort=Portfolio();
    currentPort.add([1000000],{zcb})
    ump.setCurrentPosition(currentPort);
    
    % Add an instrument to the current position. Returns a
            % resetData object so this change can be reversed easily
            resetData.currentPosition = p.currentPosition;
            resetData.payoff0 = p.payoff0;

            port = Portfolio();
            additionalQ = ones(2,1);
            additionalQ(2) = quantity;
            additionalInstruments = cell(1,2);
            additionalInstruments{1}=p.currentPosition;
            additionalInstruments{2}=instrument;
            port.add(additionalQ, additionalInstruments);
            p.currentPosition = port;         
            p.payoff0 = p.payoff0 + quantity*instrument.payoff( p.scenarios ); 
            
            
            
   % Add an instrument to the current position. Returns a
            % resetData object so this change can be reversed easily
            resetData.currentPosition = p.currentPosition;
            resetData.payoff0 = p.payoff0;

            port = Portfolio();
            port.add([p.currentPosition.map('Bond').quantity+quantity],{instrument})
            p.currentPosition = port;         
            %p.payoff0 = p.payoff0 + quantity*instrument.payoff( p.scenarios );  