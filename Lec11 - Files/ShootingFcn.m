function F = ShootingFcn (x2_0)

xx0 = [1; x2_0] ;

[~, xx] = ode45(@(t, xx) [xx(2); -0.1*xx(2)-sin(xx(1))+sin(t)], [0 20], xx0) ;

F = xx(end,2) - 0.08 ;

end
