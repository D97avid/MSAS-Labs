% Modeling and Simulation of Aerospace Systems (2021/2022)
% Assignment # 1
% Author: David Reina

%% Ex 1
clearvars; close all; clc;
format long;

% Data
f = @(x) cos(x) - x;            % Given function
a = 0; b = pi;                  % a, b parameters guessed
check = f(a)*f(b);              % Check if f(a)*f(b)<0

fprintf("f(a)*f(b) = %.4f \n",check)

% Function plot
figure()                        
xx = linspace(-5,5,1e2);
plot(xx,f(xx));
title("f(x) = cos(x) - x");
xlabel("x"); ylabel("f(x)");  
grid on

% Initialization
N = 1e5;            % Number of repetitions to obtain an average computational time value
t = zeros(N,3);     % Average computational time vector             
tol = 1e-8;         % Tolerance set (8-digit accuracy)

for i = 1:N
    % Bisection Method
    tic
    [x_BM,f_eval_BM,iter_BM] = BisecMethod(f,a,b,tol);
    t(i,1) = toc;

    % Secant Method
    tic
    [x_SM,f_eval_SM,iter_SM] = SecantMethod(f,a,b,tol);
    t(i,2) = toc;

    % Regula Falsi Method
    tic
    [x_RFM,f_eval_RFM,iter_RFM] = RegFalMethod(f,a,b,tol);
    t(i,3) = toc;
end

t_BM =  mean(t(:,1)); % Average computational time (bisection)
t_SM =  mean(t(:,2)); % Average computational time (secant)
t_RFM = mean(t(:,3)); % Average computational time (regula falsi)

RESULTS = {'Method:', 'Computational time:', 'Function evaluations:','Iterations:','Results:'; ...
           'Bisection Method',      t_BM, f_eval_BM, iter_BM, x_BM(end); ...
           'Secant Method',         t_SM, f_eval_SM, iter_SM, x_SM(end); ...
           'Regula Falsi Method',   t_RFM,f_eval_RFM,iter_RFM,x_RFM(end)};
disp(RESULTS);
format short
%% Ex 2
clearvars; close all; clc;
format long

% Data
f = @(xx) [xx(1)^2 - xx(1) - xx(2); ...     % given function
          (xx(1)^2)/16 + xx(2)^2 - 1];

% inverse of Jacobian (analytic)
invJ_an = @(xx) inv([2*xx(1) - 1, -1; 1/8 * xx(1), 2*xx(2)]);

% size of perturbation
epsilon = @(xx) max(sqrt(eps), sqrt(eps)*abs(xx));

% inverse of Jacobian (forward differences)
invJ_fd = @(xx) inv([(f(xx+[epsilon(xx(1)); 0])-f(xx))/epsilon(xx(1)), ...
                     (f(xx+[0; epsilon(xx(2))])-f(xx))/epsilon(xx(2))]);

% inverse of Jacobian (central differences)
invJ_cd = @(xx) inv([(f(xx+[epsilon(xx(1)); 0])-f(xx-[epsilon(xx(1)); 0]))/2/epsilon(xx(1)), ...
                     (f(xx+[0; epsilon(xx(2))])-f(xx-[0; epsilon(xx(2))]))/2/epsilon(xx(2))]);

% Initialization
N = 1e4;                    % Number of repetitions to obtain an average computational time value 
iter = 5;                   % number of iterations      
err_an = zeros(2,iter,2);   % vector of absolute error of analytic solution
err_fd = err_an;            % vector of absolute error of FD solution
err_cd = err_an;            % vector of absolute error of CD solution
X0 = [[-5;5],[5;5]];        % initial guesses

for j = 1:2
    x0 = X0(:,j);
    for i = 1:N
        % Analytical jacobian
        tic
        [~,F_an] = Newton(f,invJ_an,x0,iter);
        t(i,1) = toc;

        % Central differences jacobian
        tic 
        [~,F_cd] = Newton(f,invJ_cd,x0,iter);
        t(i,2) = toc;
    
        % Forward differences jacobian
        tic
        [~,F_fd] = Newton(f,invJ_fd,x0,iter);
        t(i,3) = toc;
    end

    t_an = mean(t(:,1));    % Average computational time (analytic)
    t_cd = mean(t(:,2));    % Average computational time (CD)
    t_fd = mean(t(:,3));    % Average computational time (FD)

    err_an(:,:,j) = abs(F_an); err_fd(:,:,j) = abs(F_fd); err_cd(:,:,j) = abs(F_cd);

    RESULTS = {'Method:','Error1:','Error2:','Computational time:'; ...
        'Analytical',         err_an(1,end,j),err_an(2,end,j),t_an; ...
        'Forward Differences',err_fd(1,end,j),err_fd(2,end,j),t_fd; ...
        'Central Differences',err_cd(1,end,j),err_cd(2,end,j),t_cd};
    disp(RESULTS);
end

% Find the minimum error between the three methods for each zero
[~,kA_1] = min([max(err_fd(1,end,:)),max(err_cd(1,end,:))]);
[~,kA_2] = min([max(err_fd(2,end,:)),max(err_cd(2,end,:))]);


% Plots
figure()
subplot(2,2,1)
semilogy(1:iter,err_an(1,:,1),'-o', 1:iter,err_fd(1,:,1),'-o', 1:iter,err_cd(1,:,1),'-o')
title('x_{1,A}'); xlabel('n_{iter}'); ylabel('error');
grid on
legend('Analytical solution','Forward differences','Central differences',Location='southwest')
subplot(2,2,2)
semilogy(1:iter,err_an(2,:,1),'-o', 1:iter,err_fd(2,:,1),'-o', 1:iter,err_cd(2,:,1),'-o')
title('x_{2,A}'); xlabel('n_{iter}'); ylabel('error'); 
grid on
legend('Analytical solution','Forward differences','Central differences',Location='southwest')
subplot(2,2,3)
semilogy(1:iter,err_an(1,:,2),'-o', 1:iter,err_fd(1,:,2),'-o', 1:iter,err_cd(1,:,2),'-o')
title('x_{1,B}'); xlabel('n_{iter}'); ylabel('error');
grid on
legend('Analytical solution','Forward differences','Central differences',Location='southwest')
subplot(2,2,4)
semilogy(1:iter,err_an(2,:,2),'-o', 1:iter,err_fd(2,:,2),'-o', 1:iter,err_cd(2,:,2),'-o')
title('x_{2,B}'); xlabel('n_{iter}'); ylabel('error');
grid on
legend('Analytical solution','Forward differences','Central differences',Location='southwest')
format short
%% Ex 3
clearvars; close all; clc;

% Data
odefun = @(t,x) x - t^2 + 1;                    % given function 
x_an = @(tt) tt.^2 + 2*tt + 1 - 0.5*exp(tt);    % analytic solution
x0 = 0.5;                                       % initial value        
H = [0.5, 0.2, 0.05, 0.01];                     % step size vector
tspan = [0,2];  

% initialization
N = 1e5;            % Number of repetitions to obtain an average time value
t_RK2 = {tspan(1):H(1):tspan(end),tspan(1):H(2):tspan(end),tspan(1):H(3):tspan(end),tspan(1):H(4):tspan(end)}; 
t_RK4 = t_RK2; x_RK2 = t_RK2; x_RK4 = t_RK2;
err_RK2 = t_RK2; err_RK4 = t_RK2;
t_cpu = zeros(2,4); error = t_cpu;

for i = 1:4
    for j = 1:N
        tic
        [t_RK2{i},x_RK2{i}] = RK2(odefun,tspan,x0,H(i));
        t(j,1) = toc;
        err_RK2{i} = abs(x_an(t_RK2{i})-x_RK2{i});
        tic
        [t_RK4{i},x_RK4{i}] = RK4(odefun,tspan,x0,H(i));
        t(j,2) = toc;
        err_RK4{i} = abs(x_an(t_RK4{i})-x_RK4{i});
    end

    t_cpu(1,i) = mean(t(:,1));       % RK2 average computational time for the i-th h value
    error(1,i) = err_RK2{i}(end);    % RK2 error for the i-th h value
    t_cpu(2,i) = mean(t(:,2));       % RK4 average computational time for the i-th h value
    error(2,i) = err_RK4{i}(end);    % RK4 error for the i-th h value
    
end

% Plots
figure()
subplot(2,1,1)
plot(linspace(0,2,100),x_an(linspace(0,2,100)),'--b','LineWidth',2)
hold on
plot(t_RK2{1},x_RK2{1}, t_RK2{2},x_RK2{2}, t_RK2{3},x_RK2{3}, t_RK2{4},x_RK2{4})
grid on
legend("Analytic","RK2_{h_1}","RK2_{h_2}","RK2_{h_3}","RK2_{h_4}")
title('Analytical and RK2 solutions'); xlabel('t [s]'); ylabel('x(t)');
subplot(2,1,2)
plot(linspace(0,2,1e3),x_an(linspace(0,2,1e3)),'--b','LineWidth',2)
hold on
plot(t_RK4{1},x_RK4{1}, t_RK4{2},x_RK4{2}, t_RK4{3},x_RK4{3}, t_RK4{4},x_RK4{4})
grid on
legend("Analytic","RK4_{h_1}","RK4_{h_2}","RK4_{h_3}","RK4_{h_4}")
title('Analytical and RK4 solutions'); xlabel('t [s]'); ylabel('x(t)');

figure()
plot(t_RK2{1},err_RK2{1},'--');
hold on;
plot(t_RK2{2},err_RK2{2},'--');
hold on;
plot(t_RK2{3},err_RK2{3},'--');
hold on;
plot(t_RK2{4},err_RK2{4},'--');
hold on;
grid on;
plot(t_RK4{1},err_RK4{1}, t_RK4{2},err_RK4{2}, t_RK4{3},err_RK4{3}, t_RK4{4},err_RK4{4})
legend("RK2_{h_1}","RK2_{h_2}","RK2_{h_3}","RK2_{h_4}","RK4_{h_1}","RK4_{h_2}","RK4_{h_3}","RK4_{h_4}");
title('RK absolute error'); xlabel('t [s]'); ylabel('|x_{analytic}(t) - x_{numeric}|');
ylim([0,1e-2])

figure()
loglog(t_cpu(1,:),error(1,:),'-ro');
hold on
loglog(t_cpu(2,:),error(2,:),'-bo');
grid on;
title("CPU time vs integration error");
xlabel("t_{CPU} [s]"); ylabel("\epsilon");
legend("RK2","RK4");

%% Ex 4
clearvars; close all; clc;

% Data
A = @(aa) [0 1; -1 2*cos(aa)];  % given function

% Forward operator (RK2)
F_RK2 = @(hh,aa) eye(2) + hh*(A(aa)+0.5*hh*A(aa)^2);

% Forward operator (RK4)
F_RK4 = @(hh,aa) eye(2) + hh/6 * (6*A(aa) + 3*hh*A(aa)^2 + hh^2*A(aa)^3 + hh^3/4*A(aa)^4);

options = optimset('Display','none');

% initialization
alpha = 0:pi/200:pi;
h_RK2 = zeros(length(alpha),1);
h_RK4 = h_RK2;
lambda_A = [h_RK2,h_RK2];

x0 = 2.5;  % initial guess

% Find h for alpha = [0,pi]
for i = 1:length(alpha)
    h_RK2(i) = abs(fzero(@(hh) (max(abs(eig(F_RK2(hh,alpha(i))))) - 1),x0,options));
    h_RK4(i) = abs(fzero(@(hh) (max(abs(eig(F_RK4(hh,alpha(i))))) - 1),x0,options));
    lambda_A(i,:) = eig(A(alpha(i)))';
end

% Find h for alpha = pi
disp(h_RK2(end)); disp(h_RK4(end));

% Ex. 3 data
H = [0.5, 0.2, 0.05, 0.01];
lambda_3 = 1;

% Plots
figure()
plot(real(h_RK2.*lambda_A),imag(h_RK2.*lambda_A),'b-')
xlabel("Re \{h\lambda\}"); ylabel("Im \{h\lambda\}");
title("RK2 stability domain")
grid on
axis equal

figure()
plot(real(h_RK4.*lambda_A),imag(h_RK4.*lambda_A),'r-')
hold on
plot(real(lambda_3.*H(1)),imag(lambda_3.*H(1)),'pb')
hold on
plot(real(lambda_3.*H(2)),imag(lambda_3.*H(2)),'pm')
hold on
plot(real(lambda_3.*H(3)),imag(lambda_3.*H(3)),'pg')
hold on
plot(real(lambda_3.*H(4)),imag(lambda_3.*H(4)),'pk')
legend('','','\lambda_{h = 0.5}','\lambda_{h = 0.2}','\lambda_{h = 0.05}','\lambda_{h = 0.01}')
xlabel("Re \{h\lambda\}"); ylabel("Im \{h\lambda\}");
title("RK4 stability domain")
grid on
axis equal

figure()
plot(real(h_RK2.*lambda_A),imag(h_RK2.*lambda_A),'b-')
hold on
plot(real(h_RK4.*lambda_A),imag(h_RK4.*lambda_A),'r-')
hold on
plot(real(lambda_3.*H(1)),imag(lambda_3.*H(1)),'pb')
hold on
plot(real(lambda_3.*H(2)),imag(lambda_3.*H(2)),'pm')
hold on
plot(real(lambda_3.*H(3)),imag(lambda_3.*H(3)),'pg')
hold on
plot(real(lambda_3.*H(4)),imag(lambda_3.*H(4)),'pk')
legend('RK2','','RK4','','\lambda_{h = 0.5}','\lambda_{h = 0.2}','\lambda_{h = 0.05}','\lambda_{h = 0.01}')
xlabel("Re \{h\lambda\}"); ylabel("Im \{h\lambda\}");
grid on
axis equal

save("h_RK4","h_RK4"); % save h_RK4 for exercise 7
%% Ex 5
clearvars; close all; clc;

% Data
A = @(aa) [0 1; -1 2*cos(aa)];      % given function
alpha = 0:pi/50:pi;                 % alpha vector
x0 = [1,1]';                        % initial guess
tol = [1e-3, 1e-4, 1e-5, 1e-6];     % chosen tolerances

x_an = @(t,aa) expm(A(aa)*t)*x0;    % analytic solution

% Forward operator (RK1)
F_RK1 = @(hh,aa) eye(2)+hh*A(aa);

% Forward operator (RK2)
F_RK2 = @(hh,aa) eye(2) + hh*(A(aa)+0.5*hh*A(aa)^2);

% Forward operator (RK4)
F_RK4 = @(hh,aa) eye(2) + hh/6 * (6*A(aa) + 3*hh*A(aa)^2 + hh^2*A(aa)^3 + hh^3/4*A(aa)^4);

% initialization
lambda_A = zeros(length(alpha),2);
x_an_end = zeros(2,length(alpha));
h = zeros(length(alpha),4);
h_RK1 = h; h_RK2 = h; h_RK4 = h; 
x_RK1_end = zeros(2,length(alpha)); x_RK2_end = x_RK1_end; x_RK4_end = x_RK1_end;
f_eval_RK1 = zeros(1,4); f_eval_RK2 = f_eval_RK1; f_eval_RK4 = f_eval_RK1; 

options = optimset('Display','none');

for j = 1:4
    for i = 1:length(alpha)
        x_an_end(:,i) = x_an(1,alpha(i))';
        lambda_A(i,:) = eig(A(alpha(i)))';
        h_RK1(i,j) = fzero(@(hh) norm(x_an_end(:,i)-(F_RK1(hh,alpha(i)))^(1/hh)*x0,"inf") - tol(j),5*tol(j),options);
        h_RK2(i,j) = fzero(@(hh) norm(x_an_end(:,i)-(F_RK2(hh,alpha(i)))^(1/hh)*x0,"inf") - tol(j),0.015,options);
        h_RK4(i,j) = fzero(@(hh) norm(x_an_end(:,i)-(F_RK4(hh,alpha(i)))^(1/hh)*x0,"inf") - tol(j),0.5,options);
    end
end

% func evaluation vs tol (alpha = pi)
h_RK1_pi = h_RK1(end,:);    % RK1 h for alpha = pi
h_RK2_pi = h_RK2(end,:);    % RK2 h for alpha = pi
h_RK4_pi = h_RK4(end,:);    % RK4 h for alpha = pi

for j = 1:4
    f_eval_RK1(j) = 1/h_RK1_pi(j);  % RK1 function evaluations  
    f_eval_RK2(j) = 2/h_RK2_pi(j);  % RK2 function evaluations
    f_eval_RK4(j) = 4/h_RK4_pi(j);  % RK4 function evaluations
end

% Plots
figure()
plot(real(h_RK1(:,1).*lambda_A),imag(h_RK1(:,1).*lambda_A),'m',"LineWidth",1)
hold on
plot(real(h_RK1(:,2).*lambda_A),imag(h_RK1(:,2).*lambda_A),'r',"LineWidth",1)
hold on
plot(real(h_RK1(:,3).*lambda_A),imag(h_RK1(:,3).*lambda_A),'g',"LineWidth",1)
hold on
plot(real(h_RK1(:,4).*lambda_A),imag(h_RK1(:,4).*lambda_A),'b',"LineWidth",2)
hold on
title("RK1"); xlabel("Re \{h\lambda\}"); ylabel("Im \{h\lambda\}");
legend('tol = 1e-3','','tol = 1e-4','','tol = 1e-5','','tol = 1e-6','')
grid on
axis equal

figure()
plot(real(h_RK2(:,1).*lambda_A),imag(h_RK2(:,1).*lambda_A),'m',"LineWidth",1)
hold on
plot(real(h_RK2(:,2).*lambda_A),imag(h_RK2(:,2).*lambda_A),'r',"LineWidth",1)
hold on
plot(real(h_RK2(:,3).*lambda_A),imag(h_RK2(:,3).*lambda_A),'g',"LineWidth",1)
hold on
plot(real(h_RK2(:,4).*lambda_A),imag(h_RK2(:,4).*lambda_A),'b',"LineWidth",1)
hold on
title("RK2"); xlabel("Re \{h\lambda\}"); ylabel("Im \{h\lambda\}");
legend('tol = 1e-3','','tol = 1e-4','','tol = 1e-5','','tol = 1e-6','')
grid on
axis equal

figure()
plot(real(h_RK4(:,1).*lambda_A),imag(h_RK4(:,1).*lambda_A),'m',"LineWidth",1)
hold on
plot(real(h_RK4(:,2).*lambda_A),imag(h_RK4(:,2).*lambda_A),'r',"LineWidth",1)
hold on
plot(real(h_RK4(:,3).*lambda_A),imag(h_RK4(:,3).*lambda_A),'g',"LineWidth",1)
hold on
plot(real(h_RK4(:,4).*lambda_A),imag(h_RK4(:,4).*lambda_A),'b',"LineWidth",1)
hold on
title("RK4"); xlabel("Re \{h\lambda\}"); ylabel("Im \{h\lambda\}");
legend('tol = 1e-3','','tol = 1e-4','','tol = 1e-5','','tol = 1e-6','')
grid on
axis equal

figure()
loglog(tol,f_eval_RK1(end,:),'-bo',"LineWidth",1);
hold on
loglog(tol,f_eval_RK2(end,:),'-ro',"LineWidth",1);
hold on
loglog(tol,f_eval_RK4(end,:),'-go',"LineWidth",1);
grid on
xlabel("Tolerance"); ylabel("Function evaluations")
legend('RK1','RK2','RK4')

%% Ex 6
clearvars; close all; clc;

% Data
A = @(aa) [0 1; -1 2*cos(aa)];      % given function
alpha = 0:pi/200:pi;                % alpha vector
theta = [0.1 0.3 0.4 0.7 0.9];      % vector of chosen values of theta

% Forward operator (RK2)
F_RK2 = @(hh,aa) eye(2) + hh*(A(aa)+0.5*hh*A(aa)^2);

% Operator (BI2)
B_BI2th = @(hh,aa,th) (F_RK2(-(1-th)*hh,aa))^(-1)*F_RK2((th*hh),aa);

% Initialization
h_BI2 = zeros(length(alpha),1);
lambda_A = zeros(length(alpha),2);

% h_BI2 computation
x0 = 6;         % initial guess
options = optimset('Display','none');
for j = 1:5
    for i = 1:length(alpha)
        h_BI2(i,j) = abs(fzero(@(hh) (max(abs(eig(B_BI2th(hh,alpha(i),theta(j))))) - 1),x0,options));
        lambda_A(i,:) = eig(A(alpha(i)))';
    end
end

% Plots
figure()
plot(real(h_BI2(:,1).*lambda_A),imag(h_BI2(:,1).*lambda_A),'b')
hold on
plot(real(h_BI2(:,2).*lambda_A),imag(h_BI2(:,2).*lambda_A),'r')
hold on
plot(real(h_BI2(:,3).*lambda_A),imag(h_BI2(:,3).*lambda_A),'g')
hold on
plot(real(h_BI2(:,4).*lambda_A),imag(h_BI2(:,4).*lambda_A),'m')
hold on
plot(real(h_BI2(:,5).*lambda_A),imag(h_BI2(:,5).*lambda_A),'k')
grid on
axis equal
title('Stability region (BI2)'); xlabel("Re \{h\lambda\}"); ylabel("Im \{h\lambda\}");
legend('\theta = 0.1','', '\theta = 0.3','', '\theta = 0.4','', '\theta = 0.7','', '\theta = 0.9','')

save("lambda_A","lambda_A"); save("h_BI2","h_BI2");  % save lambda_A and h_BI2 for exercise 7
%% Ex 7
clearvars; close all; clc;

% Data

B = [-180.5, 219.5; 179.5, -220.5];     % given matrix
lambda = eig(B);                        % B eigenvalues
x0 = [1,1]';                            % initial value
odefun = @(t,x) B*x;                    % given equation
an_fun = @(t) expm(B*t)*x0;             % analytic solution
h = 0.1;                                % step size
tspan = [0,5];                          % time interval
theta = 0.1;                            % given theta

% Initialization
x_an = zeros(length(tspan),2)';

% load previously calculated stability regions (BI2 and RK4)
load('lambda_A'); load('h_BI2'); load('h_RK4');

% numerical integration
[t,x_RK4] = RK4(odefun,tspan,x0,h);
[~,x_BI2] = BI2(B,tspan,x0,h,theta);

% analytic solution
for i = 1:length(t)
    x_an(:,i) = an_fun(t(i));
end

% absolute errors
errorRK4 = abs(x_an - x_RK4); 
errorBI2 = abs(x_an - x_BI2);

% Plots 
figure()
plot(real(h_RK4.*lambda_A),imag(h_RK4.*lambda_A),'r-')
hold on;
plot(real(h_BI2(:,1).*lambda_A),imag(h_BI2(:,1).*lambda_A),'b')
hold on
plot(real(lambda(1).*h),imag(lambda(1).*h),'p')
hold on
plot(real(lambda(2).*h),imag(lambda(2).*h),'p')
grid on; 
xlabel("Re \{h\lambda\}"); ylabel("Im \{h\lambda\}"); ylim([-10 10]);
title("Eigenvalue location");
legend('RK4','','BI2_{0.1}','','\lambda_1','\lambda_2')

figure()
subplot(2,1,1)
title("Numerical and analytical results comparison");
plot(t,x_an(1,:), t,x_RK4(1,:), t,x_BI2(1,:))
xlim([-1 5]); ylim([-2 2]);
xlabel("Time [s]"); ylabel("x(t)"); title("");
grid on;
legend('Analytical', 'RK4' , 'BI2_{0.1}')
subplot(2,1,2)
plot(t,x_an(2,:), t,x_RK4(2,:), t,x_BI2(2,:))
xlim([-1 5]); ylim([-1 3]);
xlabel("Time [s]"); ylabel("x(t)"); 
grid on;
legend('Analytical', 'RK4' , 'BI2_{0.1}')

figure()
subplot(2,1,1)
semilogy(t,errorBI2(1,:) , t,errorRK4(1,:))
grid on;
legend('BI2', 'RK4');
xlabel("Time [s]"); ylabel("Error");
title("x_1 absolute error")
subplot(2,1,2)
semilogy(t,errorBI2(2,:) , t,errorRK4(2,:))
grid on;
legend('BI2', 'RK4'); 
xlabel("Time [s]"); ylabel("Error");
title("x_2 absolute error")

%% Functions

function [x,f_eval,iter] = BisecMethod(fun,a,b,tol)
%
% Function performing the bisection method
%
% INPUTS:
% - fun:    function
% - a,b:    starting interval endpoints 
% - tol:    tolerance
%
% OUTPUTS:
% - x:      results vector
% - f_eval: number of function evaluations
% - iter:   number of iterations
%
funA = fun(a); funB = fun(b);
f_eval = 2; x = [];
iter = 0;
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
    iter = iter +1;
end
end

function [x,f_eval,iter] = SecantMethod(fun,a,b,tol)
%
% Function performing the secant method
%
% INPUTS:
% - fun:    function
% - a,b:    initial guesses 
% - tol:    tolerance
%
% OUTPUTS:
% - x:      results vector
% - f_eval: number of function evaluations
% - iter:   number of iterations
%
funA = fun(a); funB = fun(b);
x = [a;b]; f_eval = 2;
iter = 0;
while abs(b-a) > tol
    c = (a*funB - b*funA)/(funB-funA);
    funC = fun(c); f_eval = f_eval +1;
    a = b; funA = funB;
    b = c; funB = funC;
    x = [x;c];
    iter = iter+1;
end

end

function [x,f_eval,iter] = RegFalMethod(fun,a,b,tol)
%
% Function performing the regula falsi method
%
% INPUTS:
% - fun:    function
% - a,b:    initial guesses 
% - tol:    tolerance
%
% OUTPUTS:
% - x:      results vector
% - f_eval: number of function evaluations
% - iter:   number of iterations
%
funA = fun(a); funB = fun(b);
x = [a; b]; f_eval = 2;
iter = 0;
while abs(x(end)-x(end-1)) > tol && funA*funB < 0
    c = (a*funB - b*funA)/(funB-funA);
    funC = fun(c);
    f_eval = f_eval +1;
    
    if funA*funC < 0
        b = c;
        funB = funC;
    else
        a = c;
        funA = funC;
    end
    x = [x;c];
    iter = iter+1;
end
end

function [x,F] = Newton(fun,invJ,x0,iter)
%
% Function performing the Newton method
%
% INPUTS:
% - fun:    function
% - invJ:   inverse of Jacobian
% - x0:     initial guess 
% - iter:   maximum number of iterations
%
% OUTPUTS:
% - x:      vector of results x-coordinate
% - F:      vector of results function value (=error)
%
x = [x0, zeros(2,iter-1)];
F = x; F(:,1) = fun(x0);

for i = 1:iter-1    
    x(:,i+1) = x(:,i) - invJ(x(:,i)) * fun(x(:,i));
    F(:,i+1) = fun(x(:, i+1));
end
end

function [t,x,f_eval] = RK1(odefun,tspan,x0,h)
%
% Function performing the RK1 method
%
% INPUTS:
% - odefun:  function
% - tspan:   time interval endpoints
% - x0:      initial value 
% - h:       stepsize
%
% OUTPUTS:
% - t:       time vector
% - x:       result vector
% - f_eval:  number of function evaluations
%
t = tspan(1):h:tspan(end);
x = zeros(length(x0),length(t));
x(:,1) = x0'; 
f_eval = 0;
for i = 1:length(t)-1
    k1 = odefun(t(i),x(:,i));
    x(:,i+1) = x(:,i) + h/2 * k1;
    f_eval = f_eval+1;
end
end

function [t,x,f_eval] = RK2(odefun,tspan,x0,h)
%
% Function performing the RK2 method
% Heun's method:
% alpha = [1; 1];
% beta =  [  1,   0; 
%          1/2, 1/2];
%
% INPUTS:
% - odefun:  function
% - tspan:   time interval endpoints
% - x0:      initial value 
% - h:       stepsize
%
% OUTPUTS:
% - t:       time vector
% - x:       result vector
% - f_eval:  number of function evaluations
%

t = tspan(1):h:tspan(end);
x = zeros(length(x0),length(t));
x(:,1) = x0';
f_eval = 0;
for i = 1:length(t)-1
    k1 = odefun(t(i),x(:,i));
    k2 = odefun(t(i)+h,x(:,i)+h*k1);
    x(:,i+1) = x(:,i) + h/2 * (k1 + k2);
    f_eval = f_eval+2;
end
end

function [t,x,f_eval] = RK4(odefun,tspan,x0,h)
%
% Function performing the RK4 method
% alpha = [0.5; 0.5; 1; 1];
% beta =  [0.5,   0,   0,   0; 
%            0, 0.5,   0,   0; 
%            0,   0,   1,   0; 
%          1/6, 1/3, 1/3, 1/6];
%
% INPUTS:
% - odefun:  function
% - tspan:   time interval endpoints
% - x0:      initial value 
% - h:       stepsize
%
% OUTPUTS:
% - t:       time vector
% - x:       result vector
% - f_eval:  number of function evaluations
%
t = tspan(1):h:tspan(end);
x = zeros(length(x0),length(t));
x(:,1) = x0';
f_eval = 0;
for i = 1:length(t)-1
    k1 = odefun(t(i),x(:,i));
    k2 = odefun(t(i)+h/2,x(:,i)+h/2*k1);
    k3 = odefun(t(i)+h/2,x(:,i)+h/2*k2);
    k4 = odefun(t(i)+h,  x(:,i)+h*k3);
    x(:,i+1) = x(:,i) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
    f_eval = f_eval+4;
end
end

function [t,x] = BI2(A,tspan,x0,h,th)
%
% Function performing the BI2 method
%
% INPUTS:
% - odefun:  function
% - tspan:   time interval endpoints
% - x0:      initial value 
% - h:       stepsize
% - th:      theta
%
% OUTPUTS:
% - t:       time vector
% - x:       result vector
%
t = tspan(1):h:tspan(end);
x = zeros(length(x0),length(t));
x(:,1) = x0';

F_RK2 = @(hh) eye(2) + hh*(A+0.5*hh*A^2);
B_BI2th = @(hh,th) (F_RK2(-(1-th)*hh))^(-1)*F_RK2((th*hh));

for i = 1:length(t)-1
    x(:,i+1) = B_BI2th(h,th)*x(:,i);
end
end