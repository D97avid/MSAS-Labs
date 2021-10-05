%% Modeling and Simulation of Aerospace Systems
% IVP and BVP in Matlab
%
% F. Topputo, Politecnico di Milano

%% 1 - Calling an ode solver
f = @(t, xx) [xx(2); -0.1*xx(2)-sin(xx(1))+sin(t)] ;
xx0 = [1; -1.5] ;
[tt, xx] = ode45(f, [0 20], xx0) ;
size(tt)
size(xx)
figure(1) ; hold on ; plot(tt, xx(:,1), 'b-o', tt, xx(:,2), 'r-o') ; grid on ;
figure(2) ; plot(tt(1:end-1), diff(tt), 'k-o') ; grid on ;

%% 2 - Mass-spring-damper
clear all ;

% problem parameters
par.x0      = 0.4 ;
par.v0      = 0.2 ;
par.omega_n = 1   ; 
par.zeta    = 0.1 ; % underdamped
par.f       = @(t) 0 ; % natural motion, for now

% 2.1 - numerical integration
options = odeset ; % default options
tic ;
[tt1, xx1] = ode45(@MassDamperSpring, [0 50], [par.x0 par.v0], options, par) ; % <- see how to pass parameters to the RHS
cputime1 = toc ;
figure(1) ; hold on ; grid on ; xlabel('t') ; ylabel('x') ;
plot(tt1,xx1(:,1)) ;

% analytic solution
x_sol1 = MassDamperSpring_sol (tt1, par) ;
plot(tt1, x_sol1, 'r*') ;

% compute error on component x
err_x1 = abs(xx1(:,1) - x_sol1) ;
figure(2) ; hold on ; grid on ; xlabel('t') ; ylabel('err_x') ;
plot(tt1, err_x1, 'b') ;

% 2.2 - variable step-size in action
options = odeset(options, 'RelTol', 1e-12, 'AbsTol', 1e-12) ;
tic ;
[tt2, xx2] = ode45(@MassDamperSpring, [0 50], [par.x0 par.v0], options, par) ;
cputime2 = toc ;
x_sol2 = MassDamperSpring_sol (tt2, par) ;
err_x2 = abs(xx2(:,1) - x_sol2) ;
figure(2) ; plot(tt2, err_x2, 'r') ;

% 2.3 - Let's use another integrator (ode113)
tic
[tt3, xx3] = ode113(@MassDamperSpring, [0 50], [par.x0 par.v0], options, par) ;
cputime3 = toc ;
x_sol3 = MassDamperSpring_sol (tt3, par) ;
err_x3 = abs(xx3(:,1) - x_sol3) ;
figure(2) ; plot(tt3, err_x3, 'g') ;
%
fprintf('No. points = %d, \t MaxErr = %d, \t CPU time = %f \n', size(xx1,1), max(err_x1), cputime1) ;
fprintf('No. points = %d, \t MaxErr = %d, \t CPU time = %f \n', size(xx2,1), max(err_x2), cputime2) ;
fprintf('No. points = %d, \t MaxErr = %d, \t CPU time = %f \n', size(xx3,1), max(err_x3), cputime3) ;

% 2.4 - Let's count the number of function evaluations (use the output structure)
options = odeset ;
sol1 = ode45(@MassDamperSpring, [0 50], [par.x0 par.v0], options, par) ;
%
options = odeset(options, 'RelTol', 1e-12, 'AbsTol', 1e-12) ;
sol2 = ode45(@MassDamperSpring, [0 50], [par.x0 par.v0], options, par) ;
%
sol3 = ode113(@MassDamperSpring, [0 50], [par.x0 par.v0], options, par) ;
%
fprintf('No. points = %d, \t fcount = %d \n', size(sol1.y,2), sol1.stats.nfevals) ;
fprintf('No. points = %d, \t fcount = %d \n', size(sol2.y,2), sol2.stats.nfevals) ;
fprintf('No. points = %d, \t fcount = %d \n', size(sol3.y,2), sol3.stats.nfevals) ;

% 2.5 - forced motion, step <- passed as a parameter
par.f      = @(t) heaviside(t) ; % step
[tt4, xx4] = ode113(@MassDamperSpring, [-5 50], [par.x0 par.v0], options, par) ;
figure(1) ; plot(tt4, xx4(:,1)) ; hold on ; plot(tt4, par.f(tt4), 'r--') ;

% 2.6 - forced motion, sin
par.f      = @(t) sin(t) ; % sin
[tt5, xx5] = ode113(@MassDamperSpring, [0 50], [par.x0 par.v0], options, par) ;
figure(1) ; plot(tt5, xx5(:,1)) ; hold on ; plot(tt5, par.f(tt5), 'r--') ;

% 2.7 - forced motion, resonance
par.zeta    = 0.0 ; % no damping
[tt6, xx6] = ode113(@MassDamperSpring, [0 100], [par.x0 par.v0], options, par) ;
figure(1) ; plot(tt6, xx6(:,1), 'g') ; hold on ; plot(tt6, par.f(tt6), 'r--') ;

%% 3 - Time to pressurize an hydraulic network

clear all ;
close all ;

% pipeline
pipe.No = [1; 2; 2] ;
pipe.length = [4.5; 8.2; 5.2] ;
pipe.diameter = 0.012  ;
pipe.tickness = 0.0012 ; 

% Copute network volume
V0 = 0 ;
for i = 1 :3 
    V0 = V0 + pipe.No(i)*pi/4*pipe.diameter^2*pipe.length(i) ;
end

% 3.1 - Numerical integration of dp/dt = beta(p) * Q(t) / V0, with p(0) = 0.1 MPa
options = odeset ;
p_i = 0.1e6 ;
[tt, pp] = ode113(@ComputePdot, [0 1.6], p_i, options, V0, pipe) ;

% plots
figure(1) ; hold on ; grid on ;
xlabel('t (s)') ; ylabel('p (MPa)') ;
plot(tt, pp/1e6, 'r') ;

% 3.2 - Compute pressurizing time with 'event' option
options = odeset(options, 'event', @press_event) ;
[tt, pp, te, pe, ie] = ode113(@ComputePdot, [0 1.5], p_i, options, V0, pipe) ;
plot(te, pe/1e6, 'b*') ;
fprintf('Time to pressurize = %f (s) \n', te) ;

%% 4 - Solution of TPBVP: The shooting method
options = optimset('Display', 'Iter') ;
x2_0 = fzero(@ShootingFcn, [-1.6 -1.4], options) ;
