%% Modeling and Simulation of Aerospace Systems
% Stability of FE and BE
%
% F. Topputo, Politecnico di Milano

%% 1 - Stbility of FE
RHS = @(x, a) a*x; % xdot = f(x) = a*x
t0 = 0 ; tf = 10; h = 1;  x0 = 1 ; Nsteps = (tf-t0)/h ;
tt = linspace(t0, tf, Nsteps+1) ;
% Cases
CASES{1,1} = {-0.1}; CASES{1,2} = {'CASE A'};
CASES{2,1} = {  -1}; CASES{2,2} = {'CASE B'};
CASES{3,1} = {  -2}; CASES{3,2} = {'CASE C'};
CASES{4,1} = {  -3}; CASES{4,2} = {'CASE D'};
% plot FE stability circle
figure(10); axis equal;
xcirc = -1 + cos(0:pi/50:2*pi); ycirc = sin(0:pi/50:2*pi) ;
plot(xcirc, ycirc);
axis([-4 +4 -4 +4]) ; axis equal ;
fill(xcirc, ycirc, [0.8 0.8 0.8]) ;
% FE numerical integration
for j = 1:4
    a = cell2mat(CASES{j,1}) ;
    % Analytical solution
    xx_an = exp(a*tt) ;
    xx_nu = FE(RHS, x0, t0, tf, h, a) ;
    figure(j) ; hold on; grid on;
    plot(tt, xx_an, 'b-', tt, xx_nu, 'r-o') ;
    text(5,0.5, strcat([cellstr(CASES{j,2})])) ;
    legend('Analytical', 'Numerical') ;
    % Stability plot
    figure(10) ; hold on; grid on;
    plot(real(a*h), imag(a*h), 'k.', 'Markersize', 30) ;
    text(real(a*h), imag(a*h), strcat([cellstr(CASES{j,2})]));
end

%% 2 - Stbility of BE
% Cases
CASES{5,1} = {  -3}; CASES{5,2} = {'CASE E'};
CASES{6,1} = {  +3}; CASES{6,2} = {'CASE F'};
% plot BE stability circle
figure(10); axis equal;
th = cos(0:pi/50:2*pi); xcirc = +1 + cos(0:pi/50:2*pi); ycirc = sin(0:pi/50:2*pi) ;
plot(xcirc, ycirc); 
axis([-4 +4 -4 +4]) ; axis equal ;
fill(xcirc, ycirc, [0.8 0.8 0.8]) ;
% BE numerical integration
for j = 5:6
    a = cell2mat(CASES{j,1}) ;
    % Analytical solution
    xx_an = exp(a*tt) ;
    xx_nu = BE(RHS, x0, t0, tf, h, a) ;
    figure(j) ; hold on; grid on;
    plot(tt, xx_an, 'b-', tt, xx_nu, 'r-o') ;
    text(5,0.5, strcat([cellstr(CASES{j,2})])) ;
    legend('Analytical', 'Numerical') ;
    % Stability plot
    figure(10) ; hold on; grid on;
    plot(real(a*h), imag(a*h), 'k.', 'Markersize', 30) ;
    text(real(a*h), imag(a*h), strcat([cellstr(CASES{j,2})]));
end