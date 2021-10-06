% Modeling and Simulation of Aerospace Systems (2021/2022)
% Assignment # 1
% Author: David Reina

%% Ex 1
clearvars; close all; clc;

f = @(x) cos(x) - x;
a = 0; b = pi; tol = 1e-8;
check = f(a)*f(b);

N = 1e3;
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
f = @(xx) [xx(1)^2 - xx(1) - xx(2); xx(1)^2/16 + xx(2)^2 - 1];
invJ_an = @(xx) inv([2*xx(1) - 1, -1; 1/8 * xx(1), 2*xx(2)]);

N = 10;
x = ones(2,N);
F_an = x; F_an(:,1) = f([1;1]); F_fd = F_an;

tic
for i = 1:N-1    
    x(:,i+1) = x(:,i) - invJ_an(x(:,i)) * f(x(:,i));
    F_an(:,i+1) = f(x(:, i+1));
end
t_an = toc;

tic
for i = 1:N-1    
    x(:,i+1) = x(:,i) - invJ_an(x(:,i)) * f(x(:,i));
    F_fd(:,i+1) = f(x(:, i+1));
end
t_fd = toc;

RESULTS = {'Method:', 'Computational time:','Error:'; 'Analytical',t_an,abs(F_an(end)); 'Finite Differencies',t_fd,abs(F_fd(end))};
disp(RESULTS);
%% Ex 3
clearvars; close all; clc;
% gli input della funzione devono avere la stessa forma delle ode di matlab
x0 = 0.5;
x_analytic = @(t) t^2 + 2*t + 1 - 0.5*e^t;

%% Ex 4
clearvars; close all; clc;
A = @(aa) [0 1; -1 2*cos(aa)];
alpha = linspace(0,2*pi,1e3);
lambda = zeros(size(alpha));

for i = 1:length(alpha)
    lambda(i,:) = eig(A(alpha(i)))';
end


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
