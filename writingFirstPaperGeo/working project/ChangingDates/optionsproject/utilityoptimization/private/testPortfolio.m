function testPortfolio()

p = Portfolio();
q = zeros(3,1);
q(1) = 2;
q(2) = -0.5;
q(3) = 1;
i = cell(1,3);
i{1} = CallOption(100, NaN, NaN,Inf,Inf);
i{2} = PutOption(110, NaN, NaN, Inf,Inf);
i{3} = CallOption(200, NaN, NaN,Inf,Inf);

p.add(q, i );
assertApproxEqual( p.deltaAtInfinity(), 3.0, 0.0001 );

end

