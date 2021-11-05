% Modeling and Simulation of Aerospace Systems (2021/2022)
% Assignment # 2
% Author: David Reina

%% Exercise 1
clearvars; close all; clc;

% Data
Data.J_1 = 0.2;
Data.J_2 = 0.1;
Data.T_0 = 0.1;
Data.f = 100;
Data.tspan = 0:1/Data.f:10;

load('samples.txt');

x0 = [0; 0; 0; 0];

KB0 = [1;1];
fun = @(KB) max(norm(acceleration(x0,Data,KB(1),KB(2)) - samples(:,2:3)));
fun(KB0)
options = optimoptions('fsolve','Display','iter');
KB = fsolve(fun,KB0,options);

[acc,t,~] = acceleration(x0,Data,KB(1),KB(2));

delta = abs(acc-samples(:,2:3));

figure()
plot(samples(:,1),samples(:,2), samples(:,1),samples(:,3))
grid on
title("samples")

figure()
plot(t,acc(:,1), t,acc(:,2))
grid on
title("k = " + KB(1) + ", b = " + KB(2))

figure()
plot(t,delta(:,1), t,delta(:,2))
grid on
title("Absolute error")
xlim([t(1) t(end)])

%% Exercise 2
clearvars; close all; clc;
% Data
fluid.rho = 890;

accum.V = 0.010; % m^3
N.pi = 2.5*1e5;  % Pa
N.gamma = 1.2;
N.p0 = 21*1e5;   % Pa

accum.k = 1.12;
CV.k = 2;

B23.D = 0.018;   % m
B23.L = 2;       % m
B23.f = 0.032;

dist.k = 12;
dist.D = 0.005;  % m

piston.Dc = 0.05;  % m 
piston.Dr = 0.022; % m
piston.m = 2;
piston.x_max = 0.200; % m
piston.k = 120000; % N/m
piston.F0 = 1000; % N
piston.F = @(x) piston.F0 + piston.k*x;

B67.D = 0.018;  % m
B67.L = 15;     % m
B67.f = 0.035;

tank.p = 0.1e6;
tank.V0 = 0.01; % m^3
tank.k = 1.12;

t0 = 0;
t1 = 1;
t2 = 1.5;
tF = 3;
tspan = linspace(t0, tF,1e3);

% Starting conditions
N.V0 = accum.V * N.pi/N.p0;
fluid.V0 = accum.V - N.V0;
fluid.m = fluid.rho*fluid.V0;
fluid.p0 = N.p0;

P_A = fluid.p0*(fluid.V0/(accum.V-fluid.V0))^N.gamma;


% Branch 2-3
B23.Q = @Q_A;
B23.A = pi*B23.D^2/4;
B23.v = @(t) B23.Q(t)/B23.A;
B23.D_hyd = 4*B23.A/B23.D/pi;

% Distributor
dist.v = zeros(length(tspan),1);
dP_dist = dist.v;
for i = 1:length(tspan)
    z(i) = ramp(tspan(i),t1,t2);
    if z(i) > 0
        dist.v(i) = B23.Q(tspan(i),accum.V,fluid.V0)./((0.5*(0.5*dist.D)^2 * (2*pi*z(i) - sin(2*pi*z(i)))));
        dP_dist(i) = 0.5*dist.k*fluid.rho*dist.v(i)^2;
    end
end

% Branch 6-7
% B67.Q = Q_B; % aaaaaa
% B67.A = pi*B67.D^2/4;
% B67.v = B67.Q/B67.A;
% B67.D_hyd = 4*B67.A/B67.D/pi;

% dP_BA2 = 0.5*accum.k*fluid.rho*B23.v^2 + 0.5*CV.k*fluid.rho*B23.v^2;
% dP_B23 = 0.5*B23.f*B23.L/B23.D_hyd*fluid.rho*B23.v^2;
% dP_B67 = 0.5*B67.f*B67.L/B67.D_hyd*fluid.rho*B67.v^2;
% dP_B7T = 0.5*tank.k*fluid.rho*B67.v^2;

piston.A1 = pi*piston.Dc^2/4;
piston.A2 = piston.A1 - pi*piston.Dr^2/4; 

beta = 1.2e6; % [Pa] COPIATO DA TOPPUTO

% x = [x,v,P_4,P_5,Q_A,Q_B]
odefun = @(t,x) [x(2);
                 (x(3)*piston.A1 - x(4)*piston.A2 - piston.F(x(1)))/piston.m;
                 beta/piston.A1/x(1) * (Q_A(t,accum.V,fluid.V0)-piston.A1*x(2));
                 beta/piston.A2/(piston.x_max-x(1)) * (-x(6)-piston.A2*x(2));
                 0;
                 piston.A2 * x(2)];

x0 = [0.001; 0; 0; 0; 0; 0];
[t,x] = ode45(odefun, tspan, x0); 


% for i = 1:length(t)
%     if x(i,1) > piston.x_max
%         x(i,:) = x(k,:);
%     else
%         k = i;
%     end
% end



P = zeros(length(t),8);

% P(:,1) = P_A;
% P(:,2) = P(:,1) - dP_A2 - dP_CV;
% P(:,3) = P(:,2) - dP_B23;
% P(:,4) = P(:,3) - dP_dist;
% P(:,5) = P(:,4);
% P(:,6) = P(:,5) - 0.5*dist.k*fluid.rho*(x(:,4)./dist.A).^2;
% P(:,7) = P(:,6) - 0.5*B67.f*B67.L/B67.D_hyd*fluid.rho*(x(:,4)./B67.A).^2;
% P(:,8) = P(:,7) - 0.5*tank.k*fluid.rho*(x(:,4)./B67.A).^2;

figure()
plot(t,x(:,4))
grid on


%% Exercise 3
clearvars; close all; clc;
R_1 = 1000;
R_2 = 100;
L = 0.001;
C = 0.001;
f = 5;
Vc_0 = 1;

tspan = [0 10];
x0 = [Vc_0; 0];

% Case 1: switch closed
odefun_1 = @(t,x) [x(2);...
                   (-x(1)-x(2)*(L/R_1 + C*R_2))/(L*C+L*C*R_2/R_1)];

[t_1,x_1] = ode45(odefun_1,tspan,x0);


% Case 2: voltage source instead of switch
v = @(t) sin(2*pi*f*t)*atan(t);
v_dot = @(t) cos(2*pi*f*t)*(2*pi*f)*atan(t) + sin(2*pi*f*t)/(1+t^2);

odefun_2 = @(t,x) [x(2);...
                   (-x(1)-x(2)*(L/R_1 + C*R_2) - v(t) - v_dot(t)*L/R_1)/(L*C+L*C*R_2/R_1)];

[t_2,x_2] = ode45(odefun_2,tspan,x0);

% Plots
figure()
plot(t_1,x_1(:,1))
grid on
ylabel("V_c [V]"); xlabel("time [s]");

figure()
plot(t_2,x_2(:,1))
grid on
ylabel("V_c [V]"); xlabel("time [s]");

%% Exercise 4
clearvars; close all; clc;

Ti_0 = 20;
Ti_F = 1000;

Te = 20;
tspan = [0:60/1e3:60];

Ti = @(t) Ti_F*ramp(t,0,1);

l = [0.001 0.01 0 0.01 0.001];
k = [7.5 0.5 200 0.5 5];
c = [0 10 0 10 0];
rho = [1000 1000 1000 1000 1000];

A = 0.5^2;
R = zeros(1,5);
R(1) = l(1)./k(1)./A;
R(3) = 1/k(3);
R(5) = l(5)./k(5)./A;

% One internal node per layer
dx = l/2; 
C = A * rho .* c .* dx; 
R(2) = dx(2)./k(2)./A;
R(4) = dx(4)./k(4)./A;

x0_1 = Ti_0*ones(1,6);
dTdt_1 = @(t,x)  [((Ti(t)-x(1))/R(1) - (x(1)-x(2))/R(2))*2/C(2);
                  ((x(1)-x(2))/R(2) - (x(2)-x(3))/R(2))*1/C(2);
                  ((x(2)-x(3))/R(2) - (x(3)-x(4))/R(3))*2/C(2);
                  ((x(3)-x(4))/R(3) - (x(4)-x(5))/R(4))*2/C(4);
                  ((x(4)-x(5))/R(4) - (x(5)-x(6))/R(4))*1/C(4);
                  ((x(5)-x(6))/R(4) - (x(6)- Te )/R(5))*2/C(4)];

[t_1,T_1] = ode45(dTdt_1,tspan,x0_1);

figure()
plot(t_1,T_1);
grid on
legend("T_1","T_A","T_2","T_3","T_B","T_4")

% Two internal nodes per layer
dx = l/3; 
C = A * rho .* c .* dx; 
R(2) = dx(2)./k(2)./A;
R(4) = dx(4)./k(4)./A;

x0_2 = Ti_0*ones(1,8);

dTdt_2 = @(t,x)  [((Ti(t)-x(1))/R(1) - (x(1)-x(2))/R(2))*2/C(2);
                  ((x(1)-x(2))/R(2) - (x(2)-x(3))/R(2))*1/C(2);
                  ((x(2)-x(3))/R(2) - (x(3)-x(4))/R(2))*1/C(2);
                  ((x(3)-x(4))/R(2) - (x(4)-x(5))/R(3))*2/C(2);
                  ((x(4)-x(5))/R(3) - (x(5)-x(6))/R(4))*2/C(4);
                  ((x(5)-x(6))/R(4) - (x(6)-x(7))/R(4))*1/C(4);
                  ((x(6)-x(7))/R(4) - (x(7)-x(8))/R(4))*1/C(4);
                  ((x(7)-x(8))/R(4) - (x(8)- Te )/R(5))*2/C(4)];

[t_2,T_2] = ode45(dTdt_2,tspan,x0_2);

figure()
plot(t_2,T_2);
grid on
legend("T_1","T_A","T_B","T_2","T_3","T_C","T_D","T_4")

% figure()
% plot(t_1,Ti(t_1));
% grid on

%% Functions
function [acc,t,x] = acceleration(x0,Data,k,b)
T = @(t) Data.T_0;
odefun = @(t,x,k,b) [x(3);...
                 x(4);...
                 k*(x(2)-x(1))/Data.J_1;...
                 (-k*(x(2)-x(1)) - sign(x(4))*b*(x(4))^2 + T(t))/Data.J_2];

[t,x] = ode45(@(t,x) odefun(t,x,k,b),Data.tspan,x0);
acc = zeros(length(t),2);
for i = 1:length(t)
    acc(i,1) = k*(x(i,2)-x(i,1))/Data.J_1;
    acc(i,2) = (-k*(x(i,2)-x(i,1)) - sign(x(i,4))*b*(x(i,4))^2 + Data.T_0)/Data.J_2;
end
end

function z = ramp(t,t1,t2)
z = zeros(length(t),1);
for i = 1:length(t)
    if t(i) <= t1
        z(i) = 0;
    elseif t(i) > t1 && t(i) <= t2
        z(i) = (t(i) - t1)./(t2-t1); 
    elseif t(i) > t2
        z(i) = 1;
    end
end
end

function Q = Q_A(t,V,V0)
t1 = 1; tF = 3;
    if t <= t1
        Q = 0;
    elseif t > t1 
        Q = (V-V0)/(tF-t1); 
    end
end

