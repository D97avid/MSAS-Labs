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