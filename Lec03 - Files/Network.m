function ff = Network (yy)

% head loss coefficients
k12 = 1.122e11 ;
k23 = 1.012e11 ;
k24 = 1.234e11 ;
k35 = 1.122e11 ;
k45 = 1.320e11 ;
k56 = 1.000e11 ;

% pressure at initial, final nodes (given)
p1bar = 21e6 ;
p6bar = 20e6 ;

% states
p1 = yy(1) ;
p2 = yy(2) ;
p3 = yy(3) ;
p4 = yy(4) ;
p5 = yy(5) ;
p6 = yy(6) ;
Q12 = yy(7) ;
Q23 = yy(8) ;
Q24 = yy(9) ;
Q35 = yy(10) ;
Q45 = yy(11) ;
Q56 = yy(12) ;

% head loss in ducts (energy)
ff(1) = p1 - p2 - k12 * abs(Q12) * Q12 ;
ff(2) = p2 - p3 - k23 * abs(Q23) * Q23 ;
ff(3) = p2 - p4 - k24 * abs(Q24) * Q24 ;
ff(4) = p3 - p5 - k35 * abs(Q35) * Q35 ;
ff(5) = p4 - p5 - k45 * abs(Q45) * Q45 ;
ff(6) = p5 - p6 - k56 * abs(Q56) * Q56 ;

% internal nodes equation (continuity)
ff(7)  = Q12 - Q24 - Q23 ;
ff(8)  = Q23 - Q35 ;
ff(9)  = Q24 - Q45 ;
ff(10) = Q35 + Q45 - Q56 ;

% boundary conditions
ff(11) = p1 - p1bar ;
ff(12) = p6 - p6bar ;

end