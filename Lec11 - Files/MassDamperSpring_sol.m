function x = MassDamperSpring_sol (t, par)

omega_bar = par.omega_n*sqrt(1-(par.zeta)^2) ;
phi = atan(-(par.v0+par.zeta*par.omega_n*par.x0)/(par.x0*omega_bar)) ;
C0  = par.x0/cos(phi) ;

x = C0*exp(-par.zeta*par.omega_n*t).*cos(omega_bar*t+phi) ;

end