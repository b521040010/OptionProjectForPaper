function assertApproxEqual( a, b, e )
if ((abs(b-a)<e)~=1)
    error( 'Expected %d, actual %d',a,b);
end
end

