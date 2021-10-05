function ff = MassDamperSpring (t, xx, par)

x = xx(1) ;
v = xx(2) ;

ff(1,1) = v ;
ff(2,1) = -2*(par.zeta)*(par.omega_n)*v - (par.omega_n)^2*x + par.f(t) ;

end