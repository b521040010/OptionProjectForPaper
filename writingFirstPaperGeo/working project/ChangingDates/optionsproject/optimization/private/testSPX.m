function testSPX()
    testFindExpiryDate()
    testFindNextExpiryDate()
end

function testFindExpiryDate()

    expiry = SPX.findExpiryDate( 2016, 3 );
    assert( expiry(3) == 18 );

    expiry = SPX.findExpiryDate( 2016, 2 );
    assert( expiry(3) == 19 );

    expiry = SPX.findExpiryDate( 2016, 4 );
    assert( expiry(3) == 15 );

    expiry = SPX.findExpiryDate( 2015, 12 );
    assert( expiry(3) == 18 );

    expiry = SPX.findExpiryDate( 2015, 11 );
    assert( expiry(3) == 20 );

    expiry = SPX.findExpiryDate( 2015, 10 );
    assert( expiry(3) == 16 );

    expiry = SPX.findExpiryDate( 2015, 9 );
    assert( expiry(3) == 18 );

    expiry = SPX.findExpiryDate( 2015, 8 );
    assert( expiry(3) == 21 );

    expiry = SPX.findExpiryDate( 2015, 7 );
    assert( expiry(3) == 17 );

    expiry = SPX.findExpiryDate( 2015, 6 );
    assert( expiry(3) == 19 );

    expiry = SPX.findExpiryDate( 2015, 5 );
    assert( expiry(3) == 15 );

    expiry = SPX.findExpiryDate( 2015, 4 );
    assert( expiry(3) == 17 );


end



function testFindNextExpiryDate()

    nextExp = SPX.findNextExpiryDate([ 2016 03 20] );
    assert( nextExp(3)==15 )

    nextExp = SPX.findNextExpiryDate([ 2016 03 21] );
    assert( nextExp(3)==15 )

    nextExp = SPX.findNextExpiryDate([ 2016 04 14] );
    assert( nextExp(1)==2016 )
    assert( nextExp(2)==04 )
    assert( nextExp(3)==15 )

    nextExp = SPX.findNextExpiryDate([ 2016 06 19] );
    assert( nextExp(1)==2016 )
    assert( nextExp(2)==7 )
    assert( nextExp(3)==15 )

    % Test with quarterly expiries
    nextExp = SPX.findNextExpiryDate([ 2016 1 19], true );
    assert( nextExp(1)==2016 )
    assert( nextExp(2)==3 )
    assert( nextExp(3)==18 )

    nextExp = SPX.findNextExpiryDate([ 2016 3 18], true );
    assert( nextExp(1)==2016 )
    assert( nextExp(2)==3 )
    assert( nextExp(3)==18 )

    nextExp = SPX.findNextExpiryDate([ 2016 3 19], true );
    assert( nextExp(1)==2016 )
    assert( nextExp(2)==6 )
    assert( nextExp(3)==17 )

    nextExp = SPX.findNextExpiryDate([ 2015 8 23], true );
    assert( nextExp(1)==2015 )
    assert( nextExp(2)==9 )
    assert( nextExp(3)==18 )



end