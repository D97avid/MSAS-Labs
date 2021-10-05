function g = FindVelocity (v, parameters)

epsilon = parameters(1) ;
L   = parameters(2) ;
D   = parameters(3) ;
rho = parameters(4) ;
mu  = parameters(5) ;
dp  = parameters(6) ;

% Compute Re
Re = rho*v*D/mu ;
f  = HeadLoss(Re, epsilon) ;

% my implicit function
g = f * L/D * 1/2 * rho*v^2 - dp ;

end