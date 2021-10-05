function [value, isterminal, direction] = press_event(t, p, V0, pipe)

pn = 21e6 ; % 21 MPa 

value = pn - p ;
isterminal = 1 ;
direction = -1 ;

end