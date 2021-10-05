function xx = FE(RHS, x0, t0, tf, h, a)

Nsteps = (tf-t0)/h ;
xx = x0 ;

for k = 1 : Nsteps
    
    xx(k+1) = xx(k) + h * RHS(xx(k), a);
    
end


end