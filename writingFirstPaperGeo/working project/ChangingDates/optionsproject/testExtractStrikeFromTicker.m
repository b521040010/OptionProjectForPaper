function  testExtractStrikeFromTicker()

assert(extractStrikeFromTicker('SPX US 01/18/14 C1475 Index')==1475)
assert(extractStrikeFromTicker('SPX US 01/18/14 P1475 Index')==1475)

end

