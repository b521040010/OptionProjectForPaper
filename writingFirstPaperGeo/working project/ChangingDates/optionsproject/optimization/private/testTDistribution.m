function testTDistribution()
t = tDistribution();
t.nu = 27;
assertApproxEqual(t.integrate(1,-999999.9, 999999.9),1.0,0.0000001);
assertApproxEqual(t.integrate([0 1],-999999.9, 999999.9),0.0,0.0000001);
assertApproxEqual(t.integrate([0 0 1],-999999.9, 999999.9),1.0,0.0000001);
end

