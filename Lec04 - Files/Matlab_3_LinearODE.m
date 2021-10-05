%% Modeling and Simulation of Aerospace Systems
%
% Analytical solution of linear ODE
%
% F. Topputo, Politecnico di Milano

%% 1 - First-order model > Free response

% Define coefficients
a0 = 1; a1 = 1 ;

% Define input as function handle (zero input for now)
f = @(t) 0;

% Define generic first-order system
xdot = @(t,x) 1/a1*(f(t) - a0*x);

% Define initial condition
x0 = 1 ;

% Free response
[tt,xx] = ode45(xdot, [0 5], x0) ;

% Plot free response
plot(tt, xx) ; clear all ;

%% 2 - First-order model > Step response

% Define coefficients
a0 = 1; a1set = [0.5 1 1.5] ;

% Define input
f = @(t) heaviside(t) ;

% Define initial condition
x0 = 0 ;

hold on ; axis([0 5 0 1]) ;
for j = 1 : 3
    a1 = a1set(j);
    % Define generic first-order system
xdot = @(t,x) 1/a1*(f(t) - a0*x);
    [tt,xx] = ode45(xdot, [-1 5], x0) ;
    plot(tt, xx) ;
end

clear all ;

%% 3 - Second-order system > Overdamped

% Define mass-sping-damper
par.x0      = 0.4 ;
par.v0      = 0.2 ;
par.omega_n = 1   ; 
par.zeta    = 1.5 ; % overdamped
par.f       = @(t) 0 ; % natural motion, for now

% integrate motion
options = odeset ;
[tt1, xx1] = ode45(@MassDamperSpring, [0 10], [par.x0 par.v0], options, par) ;
plot(tt1,xx1(:,1)) ; hold on ;

%% 3 - Second-order system > Critically damped

par.zeta    = 1 ; % Critically damped

% integrate motion
[tt2, xx2] = ode45(@MassDamperSpring, [0 10], [par.x0 par.v0], options, par) ;
plot(tt2,xx2(:,1)) ;

%% 4 - Second-order system > Underdamped

par.zeta    = 0.2 ; % Underdamped

% integrate motion
[tt3, xx3] = ode45(@MassDamperSpring, [0 15], [par.x0 par.v0], options, par) ;
plot(tt3,xx3(:,1)) ;

%% 5 - Second-order systems > Step response

par.f      = @(t) 0.2* heaviside(t) ; % step
par.zeta    = 1.5 ; % overdamped

% integrate motion
[tt4, xx4] = ode45(@MassDamperSpring, [0 15], [par.x0 par.v0], options, par) ;
plot(tt4,xx4(:,1)) ; hold on ;

par.zeta    = 0.2 ; % underdamped
[tt5, xx5] = ode45(@MassDamperSpring, [0 15], [par.x0 par.v0], options, par) ;
plot(tt5,xx5(:,1)) ; hold on ;

%% 6 - Second-order systems > Harmonic input

par.f      = @(t) 0.2* sin(0.9*t) ; % input
par.zeta    = 1.8 ; % underdamped
[tt6, xx6] = ode45(@MassDamperSpring, [0 25], [par.x0 par.v0], options, par) ;
plot(tt6,xx6(:,1)) ; hold on ;

par.zeta    = 0.1 ; % underdamped
[tt7, xx7] = ode45(@MassDamperSpring, [0 25], [par.x0 par.v0], options, par) ;
plot(tt7,xx7(:,1)) ; plot(tt7, par.f(tt7)) ;

%% Normalized amplitude
zeta = [0.05, 0.1, 0.2, 0.5, 0.8] ;
omega = 0:0.001:2 ;
hold on ;
for j = 1:size(zeta,2)
    Amp = 1./sqrt([1-(omega/par.omega_n).^2].^2 + (2*zeta(j)*omega/par.omega_n).^2) ;
    plot(omega, Amp) ;
end
axis([0 2 0 8]);