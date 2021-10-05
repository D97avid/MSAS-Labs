%% Modeling and Simulation of Aerospace Systems
%
% Some notes on Matlab: Solution of nonlinear equations wih applications
%
% F. Topputo, Politecnico di Milano

%% 1 - Symbolic manipulation

% Declaration of symbolic variables
syms x y

% Definition of functions
f = cos(x) + 3*x^2
g = exp(-y^2)

% f, g are functions of pre-declared symbolic variables x, y, therefore
% they are created as symbolic functions; i.e., there is no need to declare
% them in advance

% Evaluation with numerical values or substitution of symbolic variables
subs(f,pi)
subs(f,1)

% The command subs can be used to substitute either numerical values of
% symbolic variables

% Composition and manipulation of functions
h = compose(g,f) % h(x) = g(f(x)) <- this is a function of 1 variable
k = f*g % k(x,y) = f(x) * g(y) <- this is a function of 2 varialbes, to evaluate it we have to specify the value of both variables
subs(k, [x,y], [0,1])

% Differentiation
dfdx = diff(f) % <- f(x) is a function of one variable, so there is no need to specify w.r.t. which variable f has to be differentiated 
dkdx = diff(k,x) % dk(x,y)/dx
dkdy = diff(k,y) % dk(x,y)/dy

% Integrals (indefinite)
F = int(f)
G = int(g)

% Integral (definite)
F1 = int(f,0,1)
double(F1) % to change a symbolic answer into a numerical answer; eval(F) also works
H = int(h,0,1) % this cannot be evaluated in a closed form, but we can still compute its numerical value with
double(H)

% Define polynomials, expand or simplify expressions
poly = (x-3)^5
polyexpanded = expand(poly)
polysimplified = simplify(polyexpanded)

% Find the zero of functions
solve(f) % in this case an analytical expression for the zero of f cannot be found, we have to proceed numerically (see later)
solve(g-1)
solve(polyexpanded) % finds the roots of the polynomial

% Save symbolic expressions into Matlab functions
f_handle = matlabFunction(f,x) % anonymous function (or function handle)
f_file   = matlabFunction(f,x, 'file', 'f_myfile.m')
fg_file  = matlabFunction([f; g], 'file', 'fg_myfile.m', 'vars', {x,y}) 

% Plot symbolic functions
ezplot(g)
ezplot(g,-2,2)
ezsurf(k)

%% 2 - Entering and manipulating matrices
vect1 = [1, 2, 3, 4] ; % row vector
vect2 = [1; 2; 3; 4] ; % column vector
A = rand(3,3) ;
b = ones(3,1) ;
C = zeros(3,3);

% basic formatting
format long
format compact

%% 3 - Linear algebra, x = A^{-1}*b
x1 = inv(A)*b
x2 = A\b

%% 4 - Scalar nonlinear equations

% plot the function sin(x)
xx = -5:0.01:5 ; plot(xx, sin(xx), 'b') ; hold on ; grid on ;
xlabel('$x$', 'Fontsize', 18, 'Interpreter', 'Latex') ;
ylabel('$\sin(x)$', 'Fontsize', 18, 'Interpreter', 'Latex') ;
set(gca, 'FontSize', 14) ;

% find the zeros of sin(x)
x0 = 0.5 ; plot(x0, sin(x0), 'g*') ; % first guess solution 1 (green dot)
text(x0+0.1, sin(x0), '$\leftarrow$ First guess 1', 'Fontsize', 14, 'Interpreter', 'Latex') ; % label for first guess 1
x = fzero(@(x) sin(x), x0) ; % call to fzero, function handle @(x) sin(x)
fprintf('x = %d \n', x) ; % print the solution
plot(x, sin(x), 'r*') ; % plot the solution (red dot)
text(x+0.1, sin(x), '$\leftarrow$ Solution 1', 'Fontsize', 14, 'Interpreter', 'Latex') ; % label for solution 1

% let's try with naother first guess solution
x0 = 2.5 ; plot(x0, sin(x0), 'm*') ; % first guess solution 1 (magenta dot)
text(x0+0.1, sin(x0), '$\leftarrow$ First guess solution 2', 'Fontsize', 14, 'Interpreter', 'Latex') ; % label for first guess 2
x = fzero(@(x) sin(x), x0) ; % call to fzero
fprintf('x = %d \n', x) ; % print the solution
plot(x, sin(x), 'r*') ; % plot the solution (red dot)
text(x+0.1, sin(x), '$\leftarrow$ Solution 2', 'Fontsize', 14, 'Interpreter', 'Latex') ; % label for solution 2

% passing parameters
close all ;
xx = -5:0.01:5 ; plot(xx, sin(2*xx), 'b') ; xlabel('x') ; ylabel('sin(x)') ; hold on ; grid on ;
omega = 2 ; % parameter, visible by the implicit function
x0 = 0.5 ; plot(x0, sin(omega*x0), 'g*') ;
x = fzero(@(x) sin(omega*x), x0) ;
fprintf('x = %d \n', x) ;
plot(x, sin(omega*x), 'r*') ;

% using a matlab function (when the implicit fucntion is more complex, we save it in a file)
options = optimset('TolX', 1e-14, 'Display', 'Iter') ;
x = fzero(@myfunction, x0, options, omega) ;
fprintf('x = %d \n', x) ;

% note: check help -> fzero

%% 5 - Application: Find the velocity of a fluid in a duct

% physical parameters
p1 = 21e6 ; p2 = 20.89e6 ; dp = p1-p2  ;
L = 10 ; D = 0.02 ;
rho = 930 ; nu = 30e-6 ; mu = rho*nu ;

%% 5.1 - Compute velocity (attempt #1, laminar flow, f = 64/Re)
v = dp*D^2/(32*mu*L) ;
Re = rho*v*D/mu ;
fprintf('v = %f (m/s), Re = %f (-) \n', v, Re) ;

%% 5.2 - Compute velocity (attempt #2, turbulent flow, f = 0.316/Re^0.25)
% we use symbolic manipulation
clear all ;
syms dp L D rho mu v; % define symbolic variables
Re = rho*v*D/mu ; % define Reynolds number
f = 0.316/Re^0.25 ;
fun = dp - f * L/D * 1/2 * rho * v^2 ;
v_sol = solve(fun, v) % <- solves 'fun' for v, cool feature!
% let's now define the numerical values
p1 = 21e6 ; p2 = 20.89e6 ; dp = p1-p2 ; L = 10 ; D = 0.02 ; rho = 930 ; nu = 30e-6 ; mu = rho*nu ;
% and let's evaluate the solution
v_num = eval(v_sol) ;
% consider the only real solution
v = v_num(1) ;
Re = rho*v*D/mu ; 
fprintf('v = %f (m/s), Re = %f (-) \n', v, Re) ;

%% 5.3 - Solution of the implicit function
% ok, stop guessing, let's solve the implicit function once for all
clear all ; 
p1 = 21e6 ; p2 = 20.89e6 ; dp = p1-p2 ; L = 10 ; D = 0.02 ; rho = 930 ; nu = 30e-6 ; mu = rho*nu ;
epsilon = 8e-3 ; % let's also consider the relative roughness
options = optimset('TolX', 1e-14, 'Display', 'Iter') ;
parameters = [epsilon L D rho mu dp] ; % vecotrs of parameters to pass to the implicit function
v0 = [3 4] ; % initial guess, look for solution in [3 4] m/s
v = fzero(@FindVelocity, v0, options, parameters) ;
Re = rho*v*D/mu ;
f = HeadLoss (Re, epsilon) ;
A = pi*D^2/4 ;
Q = A*v ;
K = f*L/D*1/2*rho/A^2 ;
fprintf('v = %f (m/s), Re = %f (-), f = %f \n', v, Re, f) ;
fprintf('A = %d (m^2), Q = %d (m^3/s), K = %d (kg/m^7) \n', A, Q, K) ;

% note 1: try to solve the same problem with different initial guesses for v0 (both scalars and intervals)

% note 2: try to change 'spline' into 'linear' in HeadLoss; does the result change?

%% 6 - Vector nonlinear equations
% Find the zeros of this two-dimensional function
myfun_vect = @(x) [2*x(1) -   x(2) - exp(-x(1)); ...
                    -x(1) + 2*x(2) - exp(-x(2))];
x0 = [-5; -5]; % initial guess
options=optimset('Display','iter');
[x,fval,exitflag] = fsolve(myfun_vect,x0,options) ; % call to fsolve

%% 7 - Application: Solution of hydraulic networks
clear all ;
pfg = 20.5e6*ones(6,1) ;
Qfg = 1e-3*ones(6,1) ;
yy0 = [pfg; Qfg] ;
options = optimset('TolX', 1e-10, 'Display', 'Iter') ;
yy  = fsolve(@Network, yy0, options) ;
for i = 1:6
    fprintf(strcat(['p', num2str(i), ' = %2.6f MPa \n']), yy(i)/1e6) ;
end
fprintf('Q12 = %d (m^3/s) \n', yy(7))  ;
fprintf('Q23 = %d (m^3/s) \n', yy(8))  ;
fprintf('Q24 = %d (m^3/s) \n', yy(9))  ;
fprintf('Q35 = %d (m^3/s) \n', yy(10)) ;
fprintf('Q45 = %d (m^3/s) \n', yy(11)) ;
fprintf('Q56 = %d (m^3/s) \n', yy(12)) ;