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

fid = fopen('samples.txt','rt');
samples = cell2mat(textscan(fid, '%f%f%f', 'MultipleDelimsAsOne',true, 'Delimiter','[;', 'HeaderLines',1));
fclose(fid);
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
data.accum.V = 0.01;    %[m^3]
data.accum.k = 1.12;    %[-]

data.gas.pi = 2.5e6;  %[Pa]
data.gas.gamma = 1.2;   %[-]
data.gas.p0 = 21e6;   %[Pa]
data.gas.V0 = data.accum.V * data.gas.pi/data.gas.p0;

data.fluid.rho = 890;
data.fluid.V0 = data.accum.V - data.gas.V0;
data.fluid.m = data.fluid.rho*data.fluid.V0;
data.fluid.p0 = data.gas.p0;

%accum.p = data.fluid.p0*(data.fluid.V0/(data.accum.V-data.fluid.V0))^data.gas.gamma;           %  quale dei due???
data.accum.p = data.fluid.p0*(data.fluid.V0/(data.accum.V))^data.gas.gamma;
data.CV.k = 2;              %[-]

data.B23.D = 0.018;         %[m]
data.B23.L = 2;             %[m]
data.B23.f = 0.032;         %[-]
data.B23.A = pi*data.B23.D^2/4;

data.dist.k = 12;           %[-]
data.dist.D = 0.005;        %[m]

data.piston.Dc = 0.05;      %[m] 
data.piston.Dr = 0.022;     %[m]
data.piston.m = 2;          %[kg]
data.piston.x_max = 0.200;  %[m]
data.piston.k = 120e3;      %[N/m]
data.piston.F0 = 1e3;       %[N]
data.piston.F = @(x) data.piston.F0 + data.piston.k*x;
data.piston.A1 = pi*data.piston.Dc^2/4;
data.piston.A2 = data.piston.A1 - pi*data.piston.Dr^2/4; 

data.B67.D = 0.018;  %[m]
data.B67.L = 15;     %[m]
data.B67.f = 0.035;  %[-]
data.B67.A = pi*data.B67.D^2/4;

data.tank.p = 0.1e6; %[Pa]
data.tank.V0 = 0.001; %[m^3]
data.tank.k = 1.12;  %[-]

data.t0 = 0;   %[s]
data.t1 = 1;   %[s]
data.t2 = 1.5; %[s]
data.tF = 3;   %[s]
tspan_01 = linspace(data.t0, data.t1,1e2);
tspan_12 = linspace(data.t1, data.t2,1e2);
tspan_2F = linspace(data.t2, data.tF,1e2);
tspan = linspace(data.t0, data.tF,1e2);



% alpha = @(t) 2*pi*ramp(t,t1,t2);
% dP_A4 = @(t,QA) 0.5*fluid.rho*QA*abs(QA)*((accum.k+CV.k+B23.f*B23.L/B23.D)/B23.A^2 + dist.k/(dist.D^2/8*(alpha(t)-sin(alpha(t))))^2);
% dP_T5 = @(t,QB) 0.5*fluid.rho*QB*abs(QB)*((tank.k+B67.f*B67.L/B67.D)/B67.A^2 + dist.k/(dist.D^2/8*(alpha(t)-sin(alpha(t))))^2);


x0 = [0; 0; data.accum.V; data.tank.V0];

[t,x] = ode45(@(t,x) pressurization(t,x,data), tspan, x0); 

% [t_01,x_01] = ode45(@(t,x) pressurization(t,x,data), tspan_01, x0); 
% [t_12,x_12] = ode45(@(t,x) pressurization(t,x,data), tspan_12, x_01(end,:)); 
% [t_2F,x_2F] = ode45(@(t,x) pressurization(t,x,data), tspan_2F, x_12(end,:)); 
% t = [t_01;t_12;t_2F];
% x = [x_01;x_12;x_2F];

param = zeros(length(t),12);

for i = 1:length(t)
    [xdot(:,i),param(i,:)] = pressurization(t(i),x(i,:),data);
end

% param = [P_A, P_2, P_3, P_4, P_5, P_6, P_7, P_T, Q1, Q2, alpha, F];

P_4old = data.accum.p - dP_A4(t,param(:,9),data);
P_4 = [P_4old, param(:,4)];
%%
figure()
plot(t,x(:,1))
grid on; title("x(t)")

figure()
plot(t,param(:,11))
grid on; title("alpha(t)")

figure()
plot(t,param(:,1))
hold on
plot(t,param(:,2))
hold on
plot(t,param(:,3))
hold on
plot(t,param(:,4))
hold on
plot(t,param(:,5))
hold on
plot(t,param(:,6))
hold on
plot(t,param(:,7))
hold on
plot(t,param(:,8))
grid on
legend('P_A','P_2', 'P_3', 'P_4', 'P_5', 'P_6', 'P_7', 'P_T')


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

l = [0.01 0.1 0.001 0.05 0.01];
k = [390 390 390 7.82 390];
c = [0 376.12 0 477.29 0];
rho = [8960 8960 8960 7260 8960];

A = 0.1;
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
dTdt_1 = @(t,x) [((Ti(t)-x(1))/R(1) - (x(1)-x(2))/R(2))*2/C(2);
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

dTdt_2 = @(t,x) [((Ti(t)-x(1))/R(1) - (x(1)-x(2))/R(2))*2/C(2);
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

function [xdot,paramout] = pressurization(t,xx,data)
x = xx(1);
v = xx(2);

% max x constraint:
if x <= 0
    x = 0;
end
if x <= 0 && v < 0
    x = 0;
    v = 0;
end
if x >= data.piston.x_max
    x = data.piston.x_max;
end
if x >= data.piston.x_max && v > 0
    x = data.piston.x_max;
    v = 0;
end

tol = deg2rad(90);

alpha = 2*acos(1-ramp(t,data.t1,data.t2)*2);

Q1 = data.piston.A1*v;
Q2 = -data.piston.A2*v;
P_A = data.accum.p;
P_T = data.tank.p;
P_2 = P_A - 0.5*data.fluid.rho*Q1*abs(Q1)*(data.accum.k+data.CV.k)/data.B23.A^2;
P_3 = P_2 - 0.5*data.fluid.rho*Q1*abs(Q1)*data.B23.f*data.B23.L/data.B23.D/data.B23.A^2;

if alpha > tol
    P_4 = P_A - 0.5*data.fluid.rho*Q1*abs(Q1)*((data.accum.k+data.CV.k+data.B23.f*data.B23.L/data.B23.D)/data.B23.A^2 + data.dist.k/(data.dist.D^2/8*(alpha-sin(alpha)))^2);
else 
    P_4 = P_A;
end

P_7 = P_T + 0.5*data.fluid.rho*Q2*abs(Q2)*data.tank.k/data.B67.A^2;
P_6 = P_7 + 0.5*data.fluid.rho*Q2*abs(Q2)*data.B67.f*data.B67.L/data.B67.D/data.B67.A^2;

if alpha > tol
    P_5 = P_T + 0.5*data.fluid.rho*Q2*abs(Q2)*((data.tank.k+data.B67.f*data.B67.L/data.B67.D)/data.B67.A^2 + data.dist.k/(data.dist.D^2/8*(alpha-sin(alpha)))^2);
else
    P_5 = P_T;
end

paramout = [P_A, P_2, P_3, P_4, P_5, P_6, P_7, P_T, Q1, Q2, rad2deg(alpha), data.piston.F(x)];

% ode definition:
% x = [x,v,V_acc,V_tank]
% xdot = [v;
%        ((data.accum.p-dP_A4(t,data.piston.A1*v,data))*data.piston.A1 - (data.tank.p-dP_5T(t,-data.piston.A2*v,data))*data.piston.A2 - data.piston.F(x))/data.piston.m;
%        data.piston.A1*v;
%        -data.piston.A2*v];

xdot = [v;
       (P_4*data.piston.A1 - P_5*data.piston.A2 - data.piston.F(x))/data.piston.m;
       data.piston.A1*v;
       -data.piston.A2*v];

end

function dP = dP_A4(t,Q,data)
dP = zeros(length(t),1);
alpha = 2*acos(1-ramp(t,data.t1,data.t2)*2);
%alpha = 2*pi*ramp(t,data.t1,data.t2);
tol = deg2rad(15);
for i = 1:length(t)
    if alpha(i) > tol
        dP(i) = 0.5*data.fluid.rho*Q(i)*abs(Q(i))*((data.accum.k+data.CV.k+data.B23.f*data.B23.L/data.B23.D)/data.B23.A^2 + data.dist.k/(data.dist.D^2/8*(alpha(i)-sin(alpha(i))))^2);
    else
        dP(i) = 0;
    end
end
end

function dP = dP_5T(t,Q,data)
dP = zeros(length(t),1);
alpha = 2*acos(1-ramp(t,data.t1,data.t2)*2);
tol = deg2rad(15);
for i = 1:length(t)
    if alpha(i) > tol
        dP(i) = 0.5*data.fluid.rho*Q(i)*abs(Q(i))*((data.tank.k+data.B67.f*data.B67.L/data.B67.D)/data.B67.A^2 + data.dist.k/(data.dist.D^2/8*(alpha(i)-sin(alpha(i))))^2);
    else
        dP(i) = 0;

    end
end
end