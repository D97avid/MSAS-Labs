% Modeling and Simulation of Aerospace Systems (2021/2022)
% Assignment # 1
% Author: David Reina

%% Ex 1
clearvars; close all; clc;

f = @(x) cos(x) - x;            % Given function
a = 0; b = pi;                  % a, b parameters guessed
check = f(a)*f(b); disp(check); % Check if f(a)*f(b)<0

figure()                        % Function plot
xx = linspace(-5,5,1e2);
plot(xx,f(xx),"LineWidth",4);
title("f(x) = cos(x) - x");
xlabel("x"); ylabel("y");  
grid on

N = 1e5;            % Number of repetitions to obtain an average time value
t = zeros(N,3);     % Average computational time vector             
tol = 1e-8;         % Tolerance set (8-digit accuracy)

% Bisection Method
for i = 1:N
    tic
    [x_BM,f_eval_BM,iter_BM] = BisecMethod(f,a,b,tol);
    t(i,1) = toc;
end
t_BM = mean(t(:,1));

% Secant Method
for i = 1:N
    tic
    [x_SM,f_eval_SM,iter_SM] = SecantMethod(f,a,b,tol);
    t(i,2) = toc;
end
t_SM = mean(t(:,2));

% Regula Falsi Method
for i = 1:N
    tic
    [x_RFM,f_eval_RFM,iter_RFM] = RegFalMethod(f,a,b,tol);
    t(i,3) = toc;
end
t_RFM = mean(t(:,3));

RESULTS = {'Method:', 'Computational time:', 'Function evaluations:','Iterations:','Results:'; ...
           'Bisection Method',t_BM,f_eval_BM,iter_BM,x_BM(end); ...
           'Secant Method',t_SM,f_eval_SM,iter_SM,x_SM(end); ...
           'Regula Falsi Method',t_RFM,f_eval_RFM,iter_RFM,x_RFM(end)};
disp(RESULTS);

%% Ex 2
clearvars; close all; clc;
f = @(xx) [xx(1)^2 - xx(1) - xx(2); (xx(1)^2)/16 + xx(2)^2 - 1];

invJ_an = @(xx) inv([2*xx(1) - 1, -1; 1/8 * xx(1), 2*xx(2)]);
%epsilon = 1;

epsilon = @(xx) max(sqrt(eps), sqrt(eps)*abs(xx));

invJ_fd = @(xx) inv([(f(xx+[epsilon(xx(1)); 0])-f(xx))/epsilon(xx(1)), (f(xx+[0; epsilon(xx(2))])-f(xx))/epsilon(xx(2))]);
invJ_cd = @(xx) inv([(f(xx+[epsilon(xx(1)); 0])-f(xx-[epsilon(xx(1)); 0]))/2/epsilon(xx(1)), (f(xx+[0; epsilon(xx(2))])-f(xx-[0; epsilon(xx(2))]))/2/epsilon(xx(2))]);

N = 5;
err_an = zeros(2,N,2); err_fd = err_an; err_cd = err_an;
X0 = [[-5;5],[5;5]];

for j = 1:2
    x0 = X0(:,j);
    for i = 1:1e4
        % Analytical
        tic
        [~,F_an] = Newton(f,invJ_an,x0,N);
        t(i,1) = toc;

        % Central differences
        tic
        [~,F_cd] = Newton(f,invJ_cd,x0,N);
        t(i,2) = toc;
    
        % Forward differences
        tic
        [~,F_fd] = Newton(f,invJ_fd,x0,N);
        t(i,3) = toc;
    end

    t_an = mean(t(:,1));
    t_cd = mean(t(:,2));
    t_fd = mean(t(:,3));
    err_an(:,:,j) = abs(F_an); err_fd(:,:,j) = abs(F_fd); err_cd(:,:,j) = abs(F_cd);
    RESULTS = {'Method:','Error1:','Error2:','Computational time:'; ...
        'Analytical',         err_an(1,end,j),err_an(2,end,j),t_an; ...
        'Forward Differences',err_fd(1,end,j),err_fd(2,end,j),t_fd; ...
        'Central Differences',err_cd(1,end,j),err_cd(2,end,j),t_cd};
    disp(RESULTS);
end

figure()
subplot(2,2,1)
semilogy(1:N,err_an(1,:,1),'-o', 1:N,err_fd(1,:,1),'-o', 1:N,err_cd(1,:,1),'-o')
title('x_{1,A}'); xlabel('n_{iter}'); ylabel('error');
grid on
legend('Analytical solution','Forward differences','Central differences',Location='southwest')
subplot(2,2,2)
semilogy(1:N,err_an(2,:,1),'-o', 1:N,err_fd(2,:,1),'-o', 1:N,err_cd(2,:,1),'-o')
title('x_{2,A}'); xlabel('n_{iter}'); ylabel('error'); 
grid on
legend('Analytical solution','Forward differences','Central differences',Location='southwest')
subplot(2,2,3)
semilogy(1:N,err_an(1,:,2),'-o', 1:N,err_fd(1,:,2),'-o', 1:N,err_cd(1,:,2),'-o')
title('x_{1,B}'); xlabel('n_{iter}'); ylabel('error');
grid on
legend('Analytical solution','Forward differences','Central differences',Location='southwest')
subplot(2,2,4)
semilogy(1:N,err_an(2,:,2),'-o', 1:N,err_fd(2,:,2),'-o', 1:N,err_cd(2,:,2),'-o')
title('x_{2,B}'); xlabel('n_{iter}'); ylabel('error');
grid on
legend('Analytical solution','Forward differences','Central differences',Location='southwest')

[~,kA_1] = min([err_an(1,end,1),err_fd(1,end,1),err_cd(1,end,1)]);
[~,kA_2] = min([err_an(2,end,1),err_fd(2,end,1),err_cd(2,end,1)]);
[~,kB_1] = min([err_an(1,end,2),err_fd(1,end,2),err_cd(1,end,2)]);
[~,kB_2] = min([err_an(2,end,2),err_fd(2,end,2),err_cd(2,end,2)]);

%% Ex 3
clearvars; close all; clc;
odefun = @(t,x) x - t^2 + 1; 
x_an = @(tt) tt.^2 + 2*tt + 1 - 0.5*exp(tt);
x0 = 0.5;
H = [0.5, 0.2, 0.05, 0.01];

tspan = {0,0,0,0}; x_RK2 = tspan; x_RK4 = tspan;
t_cpu = zeros(2,4); error = t_cpu;
N = 1e5;            % Number of repetitions to obtain an average time value
for i = 1:4
    tspan{i} = 0:H(i):2;
    for j = 1:N
        tic
        [~,x_RK2{i}] = RK2(odefun,tspan{i},x0,H(i));
        t(j,1) = toc;
    end
    t_cpu(1,i) = mean(t(:,1));
    error(1,i) = abs(x_an(2)-x_RK2{i}(end));
    for j = 1:N
        tic
        [~,x_RK4{i}] = RK4(odefun,tspan{i},x0,H(i));
        t(j,1) = toc;
    end
    t_cpu(2,i) = mean(t(:,1));
    error(2,i) = abs(x_an(2)-x_RK4{i}(end));
end

figure()
subplot(2,1,1)
plot(linspace(0,2,1e3),x_an(linspace(0,2,1e3)),'--b','LineWidth',2)
hold on
plot(tspan{1},x_RK2{1}, tspan{2},x_RK2{2}, tspan{3},x_RK2{3}, tspan{4},x_RK2{4})
grid on
legend("Analytic","RK2_{h_1}","RK2_{h_2}","RK2_{h_3}","RK2_{h_4}")
title('Analytical and RK2 solutions'); xlabel('t [s]'); ylabel('x(t)');

subplot(2,1,2)
plot(linspace(0,2,1e3),x_an(linspace(0,2,1e3)),'--b','LineWidth',2)
hold on
plot(tspan{1},x_RK4{1}, tspan{2},x_RK4{2}, tspan{3},x_RK4{3}, tspan{4},x_RK4{4})
grid on
legend("Analytic","RK4_{h_1}","RK4_{h_2}","RK4_{h_3}","RK4_{h_4}")
title('Analytical and RK4 solutions'); xlabel('t [s]'); ylabel('x(t)');

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
A = @(aa) [0 1; -1 2*cos(aa)];
alpha = linspace(0,2*pi,1e3);
h_RK2 = zeros(length(alpha),1);
h_RK4 = h_RK2;
lambda_A = zeros(length(alpha),2);
x0_2 = 1.2;
x0_4 = 5;

F_RK2 = @(hh,aa) eye(2) + hh*(A(aa)+0.5*hh*A(aa)^2);
F_RK4 = @(hh,aa) eye(2) + hh/6 * (6*A(aa) + 3*hh*A(aa)^2 + hh^2*A(aa)^3 + hh^3/4*A(aa)^4);

options = optimset('Display','none');
h_RK2_pi = fzero(@(hh) (max(abs(eig(F_RK2(hh,pi)))) - 1),2,options); 
h_RK4_pi = fzero(@(hh) (max(abs(eig(F_RK4(hh,pi)))) - 1),2,options);
disp(h_RK2_pi); disp(h_RK4_pi);

for i = 1:length(alpha)
    h_RK2(i) = abs(fzero(@(hh) (max(abs(eig(F_RK2(hh,alpha(i))))) - 1),x0_2,options));
    h_RK4(i) = abs(fzero(@(hh) (max(abs(eig(F_RK4(hh,alpha(i))))) - 1),x0_4,options));
    lambda_A(i,:) = eig(A(alpha(i)))';
end

figure()
plot(real(h_RK2.*lambda_A),imag(h_RK2.*lambda_A),'b-')
hold on
plot(real(h_RK4.*lambda_A),imag(h_RK4.*lambda_A),'r-')
title("Stability region")
grid on
axis equal

% Ex. 3 data
H = [0.5, 0.2, 0.05, 0.01];
lambda_3 = 1;

hold on
plot(real(lambda_3.*H),imag(lambda_3.*H),'o')
legend('RK2','','RK4')

save('h_RK4'); save('lambda_A');
%% Ex 5
clearvars; close all; clc;
A = @(aa) [0 1; -1 2*cos(aa)];

x0 = [1,1]';
alpha = [0:pi/50:2*pi];
tol = [1e-3, 1e-4, 1e-5, 1e-6];

x_an = @(t,aa) expm(A(aa)*t)*x0;
x_an_end = zeros(2,length(alpha));
for j = 1:4
    for i = 1:length(alpha)
        x_an_end(:,i) = x_an(1,alpha(i))';
        lambda_A(i,:) = eig(A(alpha(i)))';
    end
end

h = zeros(length(alpha),4);
h_RK1 = h; h_RK2 = h; h_RK4 = h; 

options = optimset('Display','none');

% RK1
x_RK1_end = zeros(2,length(alpha));
F_RK1 = @(hh,aa) eye(2)+hh*A(aa);
tic
for j = 1:4
    for i = 1:length(alpha)
        [h_RK1(i,j),~,~,out] = fzero(@(hh) norm(x_an_end(:,i)-(F_RK1(hh,alpha(i)))^(1/hh)*x0,"inf") - tol(j),5*tol(j),options);
        if h_RK1(i,j) <0
            h_RK1(i,j) = NaN;
        end
    end
end
t_RK1 = toc;

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

fprintf("RK1 done\n")

% RK2

F_RK2 = @(hh,aa) eye(2) + hh*(A(aa)+0.5*hh*A(aa)^2);
x_RK2_end = zeros(2,length(alpha));
tic
for j = 1:4
    for i = 1:length(alpha)
        [h_RK2(i,j),~,~,out] = fzero(@(hh) norm(x_an_end(:,i)-(F_RK2(hh,alpha(i)))^(1/hh)*x0,"inf") - tol(j),0.015,options);
        if h_RK2(i,j) <0
            h_RK2(i,j) = NaN;
        end
    end
end
t_RK2 = toc;

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

fprintf("RK2 done\n")

% RK4
F_RK4 = @(hh,aa) eye(2) + hh/6 * (6*A(aa) + 3*hh*A(aa)^2 + hh^2*A(aa)^3 + hh^3/4*A(aa)^4);
x_RK4_end = zeros(2,length(alpha));
tic
for j = 1:4
    for i = 1:length(alpha)
        [h_RK4(i,j),~,~,out] = fzero(@(hh) norm(x_an_end(:,i)-(F_RK4(hh,alpha(i)))^(1/hh)*x0,"inf") - tol(j),0.5,options);
        if h_RK4(i,j) <0
            h_RK4(i,j) = NaN;
        end
    end
end
t_RK4 = toc;

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

fprintf("RK4 done\n")

% func evaluation vs tol (alpha = pi -> alpha(51))
f_eval_RK1 = zeros(1,4); f_eval_RK2 = zeros(1,4); f_eval_RK4 = zeros(1,4); h_RK1_pi = h_RK1(51,:);
h_RK2_pi = h_RK2(51,:);
h_RK4_pi = h_RK4(51,:);
x0 = [1,1]';
odefun = @(t,x) A(pi)*x;

for i = 1:4
    [~,~,f_eval_RK1(i)] = RK1(odefun,0:h_RK1_pi(i):1,x0,h_RK1_pi(i));
    [~,~,f_eval_RK2(i)] = RK2(odefun,0:h_RK2_pi(i):1,x0,h_RK2_pi(i));
    [~,~,f_eval_RK4(i)] = RK4(odefun,0:h_RK4_pi(i):1,x0,h_RK4_pi(i));
end

figure()
loglog(tol,f_eval_RK1(end,:),'-bo',"LineWidth",1);
hold on
loglog(tol,f_eval_RK2(end,:),'-ro',"LineWidth",1);
hold on
loglog(tol,f_eval_RK4(end,:),'-go',"LineWidth",1);
grid on
title('Function evaluation  vs  tol')
legend('RK1','RK2','RK4')

fprintf('Elapsed time is %f seconds.',t_RK1+t_RK2+t_RK4)

%% Ex 6
clearvars; close all; clc;
theta = [0.1 0.3 0.4 0.7 0.9];
A = @(aa) [0 1; -1 2*cos(aa)];
alpha = linspace(0,2*pi,1e3);
h_BI2 = zeros(length(alpha),1);
lambda_A = zeros(length(alpha),2);
x0_2 = 6;

F_RK2 = @(hh,aa) eye(2) + hh*(A(aa)+0.5*hh*A(aa)^2);
B_RK2 = @(hh,aa) eye(2) - hh*(A(aa)-0.5*hh*A(aa)^2);
B_BI2th = @(hh,aa,th) (B_RK2((1-th)*hh,aa))^(-1)*F_RK2((th*hh),aa);

tic
options = optimset('Display','none');
for j = 1:5
    for i = 1:length(alpha)
        h_BI2(i,j) = abs(fzero(@(hh) (max(abs(eig(B_BI2th(hh,alpha(i),theta(j))))) - 1),x0_2,options));
        lambda_A(i,:) = eig(A(alpha(i)))';
    end
end
t_BI2 = toc;

disp(t_BI2); 

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
title('Stability region (BI2)')
legend('\theta = 0.1','', '\theta = 0.3','', '\theta = 0.4','', '\theta = 0.7','', '\theta = 0.9','')

% BI4
%
% h_BI4 = zeros(length(alpha),1);
% x0_4 = 14;
% F_RK4 = @(hh,aa) eye(2) + hh/6 * (6*A(aa) + 3*hh*A(aa)^2 + hh^2*A(aa)^3 + hh^3/4*A(aa)^4);
% B_RK4 = @(hh,aa) eye(2) - hh/6 * (6*A(aa) - 3*hh*A(aa)^2 + hh^2*A(aa)^3 - hh^3/4*A(aa)^4);
% B_BI4th = @(hh,aa,th) (B_RK4((1-th)*hh,aa))^(-1)*F_RK4(th*hh,aa);
% tic
% for j = 1:5
%     for i = 1:length(alpha)
%         h_BI4(i,j) = fzero(@(hh) (max(abs(eig(B_BI4th(hh,alpha(i),theta(j))))) - 1),x0_4,options);
%         lambda(i,:) = eig(A(alpha(i)))';
%         if h_BI4(i,j) < 0
%             h_BI4(i,j) = NaN;
%         end
%     end
% end
% t_BI4 = toc;
% figure()
% plot(real(h_BI4(:,1).*lambda),imag(h_BI4(:,1).*lambda),'b')
% hold on
% plot(real(h_BI4(:,2).*lambda),imag(h_BI4(:,2).*lambda),'r')
% hold on
% plot(real(h_BI4(:,3).*lambda),imag(h_BI4(:,3).*lambda),'g')
% hold on
% plot(real(h_BI4(:,4).*lambda),imag(h_BI4(:,4).*lambda),'m')
% hold on
% plot(real(h_BI4(:,5).*lambda),imag(h_BI4(:,5).*lambda),'k')
% grid on
% axis equal
% title('Stability region (BI4)')
% legend('\theta = 0.1','', '\theta = 0.3','', '\theta = 0.4','', '\theta = 0.7','', '\theta = 0.9','')
%disp(t_BI4);

save('h_BI2');
%% Ex 7
clearvars; close all; clc;
% NB B is a STIFF MATRIX
B = [-180.5, 219.5; 179.5, -220.5];
x0 = [1,1]';
h = 0.1; tspan = 0:h:5;
theta = 0.1;

% analytic solution
x_an = zeros(length(tspan),2)';
an_fun = @(t) expm(B*t)*x0;
odefun = @(t,x) B*x;
for i = 1:length(tspan)
    x_an(:,i) = an_fun(tspan(i));
end

% numerical integration
[~,x_RK4] = RK4(odefun,tspan,x0,h);
[~,x_BI2] = BI2(B,tspan,x0,h,theta);

% h-lambda plane
lambda = eig(B);
load('lambda_A'); load('h_BI2'); save('h_RK4');

figure()
plot(real(h_RK4.*lambda_A),imag(h_RK4.*lambda_A),'r-')
hold on;
plot(real(h_BI2(:,1).*lambda_A),imag(h_BI2(:,1).*lambda_A),'b')
hold on
plot(real(lambda.*h),imag(lambda.*h),'o')
grid on; axis equal;
legend('RK4','','BI4','')

figure()
subplot(2,1,1)
plot(tspan,x_an(1,:), tspan,x_BI2(1,:))
grid on
axis equal
legend('Analytical', 'BI2')
subplot(2,1,2)
plot(tspan,x_an(2,:), tspan,x_BI2(2,:))
grid on
axis equal
legend('Analytical', 'BI2')

figure()
subplot(2,1,1)
plot(tspan,x_an(1,:), tspan,x_RK4(1,:))
xlim([0 5]); ylim([-2 2]);
grid on; axis equal;
legend('Analytical', 'RK4')

subplot(2,1,2)
plot(tspan,x_an(2,:), tspan,x_RK4(2,:))
xlim([0 5]); ylim([-2 2]);
grid on; axis equal;
legend('Analytical', 'RK4')

%% Functions
function [x,f_eval,iter] = BisecMethod(fun,a,b,tol)
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
x = [a;b]; f_eval = 0;
iter = 0;
while abs(b-a) > tol
    funA = fun(a); funB = fun(b);
    f_eval = f_eval +2;
    c = (a*funB - b*funA)/(funB-funA);
    a = b; b = c;
    x = [x;c];
    iter = iter+1;
end

end

function [x,f_eval,iter] = RegFalMethod(fun,a,b,tol)
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

function [x,F] = Newton(f,invJ,x0,N)
x = [x0, zeros(2,N-1)];
F = x; F(:,1) = f(x0);

for i = 1:N-1    
        x(:,i+1) = x(:,i) - invJ(x(:,i)) * f(x(:,i));
        F(:,i+1) = f(x(:, i+1));
end
end

function [tspan,x,f_eval] = RK1(odefun,tspan,x0,h)
x = zeros(length(x0),length(tspan));
x(:,1) = x0'; 
f_eval = 0;
for i = 1:length(tspan)-1
    k1 = odefun(tspan(i),x(:,i));
    x(:,i+1) = x(:,i) + h/2 * k1;
    f_eval = f_eval+1;
end
end

function [tspan,x,f_eval] = RK2(odefun,tspan,x0,h)
% a = [1; 1];
% b = [1, 0; 1/2, 1/2];
x = zeros(length(x0),length(tspan));
x(:,1) = x0';
f_eval = 0;
for i = 1:length(tspan)-1
    k1 = odefun(tspan(i),x(:,i));
    k2 = odefun(tspan(i)+h,x(:,i)+h*k1);
    x(:,i+1) = x(:,i) + h/2 * (k1 + k2);
    f_eval = f_eval+2;
end
end

function [tspan,x,f_eval] = RK4(odefun,tspan,x0,h)
% a = [0.5; 0.5; 1; 1];
% b = [0.5, 0, 0, 0; 0, 0.5, 0, 0; 0, 0, 1, 0; 1/6, 1/3, 1/3, 1/6];
x = zeros(length(x0),length(tspan));
x(:,1) = x0';
f_eval = 0;
for i = 1:length(tspan)-1
    k1 = odefun(tspan(i),x(:,i));
    k2 = odefun(tspan(i)+h/2,x(:,i)+h/2*k1);
    k3 = odefun(tspan(i)+h/2,x(:,i)+h/2*k2);
    k4 = odefun(tspan(i)+h,  x(:,i)+h*k3);
    x(:,i+1) = x(:,i) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
    f_eval = f_eval+4;
end
end

function [tspan,x] = BI2(A,tspan,x0,h,th)

F_RK2 = @(hh) eye(2) + hh*(A+0.5*hh*A^2);
B_RK2 = @(hh) eye(2) - hh*(A-0.5*hh*A^2);
B_BI2th = @(hh,th) (B_RK2((1-th)*hh))^(-1)*F_RK2((th*hh));

x = zeros(length(x0),length(tspan));
x(:,1) = x0';

for i = 1:length(tspan)-1
    x(:,i+1) = B_BI2th(h,th)*x(:,i);
end
end