function pdot = ComputePdot (t, p, V0, pipe)

pn = 21e6 ; % 21 MPa 
Wn = 10e3 ; % 10 kW
tr = 5    ; %  5 s
E  = 70e9 ; % 70 GPa
gas_conc = 0.002 ; % 0.2%, gas concentration

% power
W = Wn/tr*t ;
% flow rate
Q = W/pn ;

% bulk modulus (fluid)
beta0 = 1200e6 ; % 1200 MPa
p0 = 0.1e6 ;
T0 = 40    ;
PI = 260e6 ;
OM = 192   ;
T  = 40    ;
beta_l = beta0*10^((p-p0)/PI - (T-T0)/OM) ;

% bulk modulus (gas)
beta_g = p ;

% bulk modulus (pipes)
beta_c = E*pipe.tickness/pipe.diameter ;

% effective bulk modulus
inv_beta_e = 1/beta_l + gas_conc * 1/beta_g + 1/beta_c ;
beta_e = 1/inv_beta_e ;

% compute dp/dt
pdot = beta_e * Q / V0 ;

end