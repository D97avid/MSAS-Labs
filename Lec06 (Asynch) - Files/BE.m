function xx = BE(RHS, x0, t0, tf, h, a)

Nsteps = (tf-t0)/h ;
xx = x0 ;

for k = 1 : Nsteps
    
    xx(k+1) = fzero(@(xkp1) xkp1 - xx(k) - h * RHS(xkp1, a) , xx(k)) ;
    
end


end