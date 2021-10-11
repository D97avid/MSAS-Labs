% Modeling and Simulation of Aerospace Systems (2021/2022)
% Assignment # 1
% Author: David Reina

%% Ex 1
clearvars; close all; clc;

f = @(x) cos(x) - x;
a = 0; b = pi; tol = 1e-8;
check = f(a)*f(b);

figure()
xxx = linspace(-5,5,1e2);
plot(xxx,f(xxx))
grid on

N = 1e5;
t = ones(N,1);

% Bisection Method
for i = 1:N
    tic
    [x_BM,f_eval_BM] = BisecMethod(f,a,b,tol);
    t(i) = toc;
end
t_BM = mean(t);

% Secant Method
for i = 1:N
    tic
    [x_SM,f_eval_SM] = SecantMethod(f,a,b,tol);
    t(i) = toc;
end
t_SM = mean(t);

% Regula Falsi Method
for i = 1:N
    tic
    [x_RFM,f_eval_RFM] = RegFalMethod(f,a,b,tol);
    t(i) = toc;
end
t_RFM = mean(t);

RESULTS = {'Method:', 'Computational time:', 'Function evaluations:','Results:'; 'Bisection Method',t_BM,f_eval_BM,x_BM(end); 'Secant Method',t_SM,f_eval_SM,x_SM(end); 'Regula Falsi Method',t_RFM,f_eval_RFM,x_RFM(end)};
disp(RESULTS);


%% Ex 2
clearvars; close all; clc;
% we don't have to impose an error tolerance but stop the iterations after
% a number of steps and compare the accuracy. we decide a priori depending
% on the newton method 
f = @(xx) [xx(1)^2 - xx(1) - xx(2); (xx(1)^2)/16 + xx(2)^2 - 1];

figure()
xxx = linspace(-3,3,1e2);
plot(xxx,(xxx.^2 - xxx) , xxx,(-xxx.^2 ./16 + 1))
grid on

invJ_an = @(xx) inv([2*xx(1) - 1, -1; 1/8 * xx(1), 2*xx(2)]);
epsilon = 1e-3;
% eps = @(xx) max(sqrt(epsilon),max(sqrt(epsilon)*abs(xx)));
% % H = diag([eps,eps]);
% invJ_fd = @(xx) inv(eps(xx)^-1 .* [f(xx+[eps(xx); 0])-f(xx), f(xx+[0; eps(xx)])-f(xx)]);
% invJ_cd = @(xx) inv((2*eps(xx))^-1 .* [f(xx+[eps(xx); 0])-f(xx-[eps(xx); 0]), f(xx+[0; eps(xx)])-f(xx-[0; eps(xx)])]);

invJ_fd = @(xx) inv(epsilon^-1 .* [f(xx+[epsilon; 0])-f(xx), f(xx+[0; epsilon])-f(xx)]);
invJ_cd = @(xx) inv((2*epsilon)^-1 .* [f(xx+[epsilon; 0])-f(xx-[epsilon; 0]), f(xx+[0; epsilon])-f(xx-[0; epsilon])]);

N = 10;
X0 = [[-5;5],[5;5]];
for j = 1:2
    x0 = X0(:,j);
    x = [x0, zeros(2,N-1)]; x_an = x; x_fd = x; x_cd = x;
    F_an = x; F_an(:,1) = f([1;1]); F_fd = F_an; F_cd = F_an;
    
    % Analytical
    tic
    for i = 1:N-1    
        x_an(:,i+1) = x_an(:,i) - invJ_an(x_an(:,i)) * f(x_an(:,i));
        F_an(:,i+1) = f(x_an(:, i+1));
    end
    t_an = toc;
    
    % Forward differences
    tic
    for i = 1:N-1
        x_fd(:,i+1) = x_fd(:,i) - invJ_fd(x_fd(:,i)) * f(x_fd(:,i));
        F_fd(:,i+1) = f(x_fd(:, i+1));
    end
    t_fd = toc;
    
    % Central differences
    tic
    for i = 1:N-1
        x_cd(:,i+1) = x_cd(:,i) - invJ_cd(x_cd(:,i)) * f(x_cd(:,i));
        F_cd(:,i+1) = f(x_cd(:, i+1));
    end
    t_cd = toc;
    
    figure()
    subplot(2,1,1)
    plot(1:N,F_an(1,:),'-o', 1:N,F_fd(1,:),'-o', 1:N,F_cd(1,:),'-o')
    grid on
    legend('Analytical solution','Forward differences','Central differences')
    subplot(2,1,2)
    plot(1:N,F_an(2,:),'-o', 1:N,F_fd(2,:),'-o', 1:N,F_cd(2,:),'-o')
    grid on
    legend('Analytical solution','Forward differences','Central differences')
    
    RESULTS = {'Method:', 'Computational time:','Error1:','Error2:','x1:','x2:'; 'Analytical',t_an,abs(F_an(1,end)),abs(F_an(2,end)),x_an(1,end),x_an(2,end); 'Forward Differences',t_fd,abs(F_fd(1,end)),abs(F_fd(2,end)),x_fd(1,end),x_fd(2,end); 'Central Differences',t_cd,abs(F_cd(1,end)),abs(F_cd(2,end)),x_cd(1,end),x_cd(2,end)};
    disp(RESULTS);
end

%% Ex 3
clearvars; close all; clc;
odefun = @(t,x) x - t^2 + 1; 
x_an = @(tt) tt.^2 + 2*tt + 1 - 0.5*exp(tt);
x0 = 0.5;
tspan = linspace(0,2,1e1)';
x_RK2 = zeros(length(tspan),4);
x_RK4 = x_RK2;
H = [0.5, 0.2, 0.05, 0.01];

for i = 1:4
    [~,x_RK2(:,i)] = RK2(odefun,tspan,x0,H(i));
    [~,x_RK4(:,i)] = RK4(odefun,tspan,x0,H(i));
end

figure()
plot(tspan,x_an(tspan),'LineWidth',3)
hold on
plot(tspan,x_RK2(:,1), tspan,x_RK2(:,2), tspan,x_RK2(:,3), tspan,x_RK2(:,4))
hold on
plot(tspan,x_RK4(:,1), tspan,x_RK4(:,2), tspan,x_RK4(:,3), tspan,x_RK4(:,4))
grid on
legend("Analytic","RK2_{h_1}","RK2_{h_2}","RK2_{h_3}","RK2_{h_4}",...
    "RK4_{h_1}","RK4_{h_2}","RK4_{h_3}","RK4_{h_4}")

%% Ex 4
clearvars; close all; clc;
A = @(aa) [0 1; -1 2*cos(aa)];
alpha = linspace(0,pi,1e3);
h_RK2 = zeros(length(alpha),1);
h_RK4 = h_RK2;
lambda = zeros(length(alpha),2);
x0 = 10;

coeff_a = [1;1]; coeff_b = [1, 0 ; 0.5, 0.5]; % Heun's method
F_RK2 = @(hh,aa) eye(2) + coeff_a(2)*hh*(coeff_b(2,1)*A(aa)+coeff_b(2,2)*A(aa)+coeff_b(1,1)*coeff_b(2,2)*hh*A(aa)^2);
F_RK4 = @(hh,aa) eye(2) + hh/6 * (6*A(aa) + 3*hh*A(aa)^2 + hh^2*A(aa)^3 + hh^3/4*A(aa)^4);

options = optimoptions('fsolve','Display','none');
h_RK2_pi = fsolve(@(hh) (max(abs(eig(F_RK2(hh,pi)))) - 1),x0,options);
h_RK4_pi = fsolve(@(hh) (max(abs(eig(F_RK4(hh,pi)))) - 1),x0,options);

for i = 1:length(alpha)
    h_RK2(i) = fsolve(@(hh) (max(abs(eig(F_RK2(hh,alpha(i))))) - 1),x0,options);
    h_RK4(i) = fsolve(@(hh) (max(abs(eig(F_RK4(hh,alpha(i))))) - 1),x0,options);
    lambda(i,:) = eig(A(alpha(i)))';
end

figure()
plot(real(h_RK2.*lambda),imag(h_RK2.*lambda))
hold on
plot(real(h_RK4.*lambda),imag(h_RK4.*lambda))
title("Stability region")
grid on
axis equal

% Ex. 3 data
H = [0.5, 0.2, 0.05, 0.01];
lambda_3 = 1;

hold on
plot(real(lambda_3.*H),imag(lambda_3.*H),'o')


%% Ex 5
clearvars; close all; clc;


%% Ex 6
clearvars; close all; clc;


%% Ex 7
clearvars; close all; clc;
% NB B is a STIFF MATRIX
% theta = h = 0.1
B = [-180.5, 219.5; 179.5, -220.5];
x0 = [1,1]';


%% Functions
function [x,f_eval] = BisecMethod(fun,a,b,tol)
x = [];
funA = fun(a);
funB = fun(b);
f_eval = 2;

while abs(b-a) > tol && funA*funB < 0
    x0 = (a+b)/2;
    x = [x; x0];
    fun0 = fun(x0);   
    f_eval = f_eval +1;

    if funA*fun0 < 0
        b = x0;
        funB = fun0;
    else
        a = x0;
        funA = fun0;
    end
end
end

function [x,f_eval] = SecantMethod(fun,a,b,tol)
x = [];
f_eval = 0;

while abs(b-a) > tol
    funA = fun(a);
    funB = fun(b);
    f_eval = f_eval +2;

    c = (a*funB - b*funA)/(funB-funA);
    a = b;
    b = c;
    x = [x;c];
end

end

function [x,f_eval] = RegFalMethod(fun,a,b,tol)

funA = fun(a);
funB = fun(b);
c = (a*funB - b*funA)/(funB-funA);
x = [c];
f_eval = 2;


while abs(b-a) > tol && funA*funB < 0
    c = (a*funB - b*funA)/(funB-funA);
    funC = fun(c);
    f_eval = f_eval +1;
    
    if funA*funC > 0
        a = c;
        funA = funC;
    else
        b = c;
        funB = funC;
    end
    x = [x;c];
end
end

function [tspan,x] = RK2(odefun,tspan,x0,h)
% a = [1; 1];
% b = [1, 0; 1/2, 1/2];
x = zeros(size(tspan));
x(1) = x0;

for i = 1:length(tspan)-1
    k1 = odefun(tspan(i),x(i));
    k2 = odefun(tspan(i)+h,x(i)+h*k1);
    x(i+1) = x(i) + h/2 * (k1 + k2);
end
end


function [tspan,x] = RK4(odefun,tspan,x0,h)
% a = [0.5; 0.5; 1; 1];
% b = [0.5, 0, 0, 0; 0, 0.5, 0, 0; 0, 0, 1, 0; 1/6, 1/3, 1/3, 1/6];
x = zeros(size(tspan));
x(1) = x0;

for i = 1:length(tspan)-1
    k1 = odefun(tspan(i),x(i));
    k2 = odefun(tspan(i)+h/2,x(i)+h/2*k1);
    k3 = odefun(tspan(i)+h/2,x(i)+h/2*k2);
    k4 = odefun(tspan(i)+h,  x(i)+h*k3);
    x(i+1) = x(i) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end
end