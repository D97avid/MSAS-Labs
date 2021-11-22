% Modeling and Simulation of Aerospace Systems (2021/2022)
% Assignment # 2
% Author: David Reina

%% Exercise 1
clearvars; close all; clc;

% Data
Data.J_1 = 0.2;                 %[kg m^2] Disk 1 moment of inertia 
Data.J_2 = 0.1;                 %[kg m^2] Disk 2 moment of inertia
Data.T_0 = 0.1;                 %[N m]    Motor torque
Data.f = 100;                   %[Hz]     Sampling frequency
Data.tspan = 0:1/Data.f:10;     %[s]      Time vetor
x0 = [0; 0; 0; 0];              % Initial conditions vector
k_guess = 1; b_guess = 1;       % Guessed k and b 

% Load file samples.txt:
fid = fopen('samples.txt','rt');
samples = cell2mat(textscan(fid, '%f%f%f', 'MultipleDelimsAsOne',true, 'Delimiter','[;', 'HeaderLines',1));
fclose(fid);

% Find system response for guessed k and b:
[acc_guess,t_guess,~] = acceleration(x0,Data,k_guess,b_guess);

% Find k and b such that parametric errors are minimized:
fun = @(KB) max(norm(acceleration(x0,Data,KB(1),KB(2)) - samples(:,2:3)));  % Function to minimize (maximum absolute error)

KB0 = [k_guess;b_guess];                                                    % Initial guess vector for fmincon
options = optimoptions("fminunc","Display","off");
[KB] = fminunc(fun,KB0,options);

k_opt = KB(1);
b_opt = KB(2);

fprintf("k_opt = %.4f\n",k_opt); fprintf("b_opt = %.4f\n",b_opt);

% Find system response for guessed k_opt and b_opt:
[acc,t,~] = acceleration(x0,Data,k_opt,b_opt);
eps = abs(acc-samples(:,2:3));                                              % Absolute error 

% Plots:
figure()
plot(t_guess,acc_guess(:,1),"LineWidth",2,"Color",[0.3010 0.7450 0.9330])
hold on
plot(t_guess,acc_guess(:,2),"LineWidth",2,"Color","#FFE42B")
grid on; legend('$\dot\theta_1$','$\dot\theta_2$','Interpreter','latex')
title('System response for $k_{guess}$ and $b_{guess}$','Interpreter','latex'); 
xlabel('$Time \,\, [s]$','Interpreter','latex'); ylabel('$\dot\theta \,\, [rad/s^2]$','Interpreter','latex');

figure()
subplot(1,2,1)
plot(samples(2:end,1),samples(2:end,2),"LineWidth",2,"LineStyle","-","Color",[0.3010 0.7450 0.9330])
hold on
plot(t(2:end),acc(2:end,1),"LineWidth",2,"LineStyle",":","Color",[0 0.4470 0.7410])
hold on
plot(samples(2:end,1),samples(2:end,3),"LineWidth",2,"LineStyle","-","Color","#FFE42B")
hold on
plot(t(2:end),acc(2:end,2),"LineWidth",2,"LineStyle",":","Color",[0.9290 0.6940 0.1250])
grid on; 
legend('$\dot\theta_{1,sampled}$','$\dot\theta_1$','$\dot\theta_{2,sampled}$','$\dot\theta_2$','Interpreter','latex')
title('Accelerations','Interpreter','latex');
xlabel('$Time \,\, [s]$','Interpreter','latex'); ylabel('$\dot\theta \,\, [rad/s^2]$','Interpreter','latex');

subplot(1,2,2)
semilogy(t(2:end),eps(2:end,1),"LineWidth",2,"Color",[0.3010 0.7450 0.9330])
hold on
semilogy(t(2:end),eps(2:end,2),"LineWidth",2,"Color","#FFE42B")
grid on; legend('$\epsilon_1$','$\epsilon_2$','Interpreter','latex')
title('Absolute error','Interpreter','latex');
xlabel('$Time \,\, [s]$','Interpreter','latex'); ylabel('$\epsilon \,\, [-]$','Interpreter','latex');

%% Exercise 2
clearvars; close all; clc;

% Data
data.accum.V = 0.01;                                      %[m^3] Total volume (accumulator)
data.accum.k = 1.12;                                      %[-] Pressure drop coefficient (accumulator)

data.gas.pi = 2.5e6;                                      %[Pa] Initial pressure (nitrogen)
data.gas.gamma = 1.2;                                     %[-] Adiabatic exponent (nitrogen)
data.gas.p0 = 21e6;                                       %[Pa] Pressure at t0 (nitrogen)
data.gas.V0 = data.accum.V * data.gas.pi/data.gas.p0;     %[m^3] Volume at t0 (nitrogen)

data.fluid.rho = 890;                                     %[kg/m^3] Density (skydrol)
data.fluid.V0 = data.accum.V - data.gas.V0;               %[m^3] Volume at t0 (skydrol)
data.fluid.p0 = data.gas.p0;                              %[Pa] Pressure at t0 (skydrol)

data.CV.k = 2;                                            %[-] Pressure drop coefficient (check valve)
                              
data.B23.D = 0.018;                                       %[m] Diameter (branch 2-3)
data.B23.L = 2;                                           %[m] Length (branch 2-3)
data.B23.f = 0.032;                                       %[-] Friction factor (branch 2-3)
data.B23.A = pi*data.B23.D^2/4;                           %[m^2] Cross section (branch 2-3)

data.dist.k = 12;                                         %[-] Pressure drop coefficient (distributor)
data.dist.D = 0.005;                                      %[m] Cross section (distributor)
                              
data.piston.Dc = 0.05;                                    %[m] Piston head diameter
data.piston.Dr = 0.022;                                   %[m] Piston rod diameter
data.piston.m = 2;                                        %[kg] Piston mass
data.piston.x_max = 0.200;                                %[m] Piston maximum stroke
data.piston.k = 120e3;                                    %[N/m] Piston load elastic coefficient
data.piston.F0 = 1e3;                                     %[N] Piston pre-load
data.piston.F = @(x) data.piston.F0 + data.piston.k*x;    %[N] Piston load function
data.piston.A1 = pi*data.piston.Dc^2/4;                   %[m^2] Piston head area
data.piston.A2 = data.piston.A1 - pi*data.piston.Dr^2/4;  %[m^2] Piston rod area

data.B67.D = 0.018;                                       %[m] Diameter (branch 6-7)
data.B67.L = 15;                                          %[m] Length (branch 6-7)
data.B67.f = 0.035;                                       %[-] Friction factor (branch 6-7)
data.B67.A = pi*data.B67.D^2/4;                           %[m^2] Cross section (branch 6-7)

data.tank.p0 = 0.1e6;                                     %[Pa] Pressure at t0 (tank)
data.tank.V0 = 0.001;                                     %[m^3] Volume at t0 (nitrogen)
data.tank.k = 1.12;                                       %[-] Pressure drop coefficient (tank)
                                    
data.t0 = 0;                                              %[s] Simulation start
data.t1 = 1;                                              %[s] Start distributor opening
data.t2 = 1.5;                                            %[s] End distributor opening
data.tF = 3;                                              %[s] Simulation end

% Integration until event (piston position = maximum stroke)
% State vector: x = [x, v, V_accumulator, V_tank]
tspan = linspace(data.t0, data.tF,1e3)';                 % Time vector
x0 = [0; 0; data.fluid.V0; data.tank.V0];                 % Initial conditions
options_BE = odeset('event', @press_event); 
[t_BE,x_BE,t_event,x_event,~] = ode23s(@pressurization_BE, tspan, x0, options_BE,data); 
param_BE = zeros(length(t_BE),13);
for i = 1:length(t_BE)
    [~,param_BE(i,:)] = pressurization_BE(t_BE(i),x_BE(i,:),data);
end
fprintf("t_event = %f\n",t_event)

% After the event, the parameter don't change except the piston velocity
% which goes to zero:
t = ones(2*length(t_BE),1);
t(1:length(t_BE),1) = t_BE;
t(length(t_BE)+1:end,1) = linspace(t_BE(end),data.tF,length(t)-length(t_BE));
x = ones(length(t),4);
x(1:length(t_BE),:) = x_BE;
x_PE = [1 0 1 1] .* x_event;
x(length(t_BE)+1:end,:) = x(length(t_BE)+1:end,:) .*x_PE;

param = ones(length(t),width(param_BE));
param(1:length(t_BE),:) = param_BE;
param(length(t_BE)+1:end,:) = param(length(t_BE)+1:end,:).*param_BE(end,:);

% Plots:
figure()
subplot(2,2,1)
plot(t,x(:,1),'LineWidth',2)
hold on
line(xlim, [0,0], 'Color', 'k', 'LineWidth', 0.5);
plot(t,x(:,1),'LineWidth',2,'Color',[0, 0.4470, 0.7410])
hold on
plot(t_event,x_event(1),'MarkerSize',10,'Marker','pentagram','MarkerFaceColor','m','MarkerEdgeColor','k',"Color","white")
grid on; title('x(t)','Interpreter','latex')
xlabel('$Time \,\, [s]$','Interpreter','latex'); ylabel('$x \,\, [m]$','Interpreter','latex');
ylim([-0.05 0.25])
legend("","","","Event",'Interpreter','latex');
subplot(2,2,2)
plot(t,x(:,2),'LineWidth',2)
grid on; title('v(t)','Interpreter','latex')
xlabel('$Time \,\, [s]$','Interpreter','latex'); ylabel('$v \,\, [m/s^2]$','Interpreter','latex');
line(xlim, [0,0], 'Color', 'k', 'LineWidth', 0.5);
hold on
plot(t,x(:,2),'LineWidth',2,'Color',[0, 0.4470, 0.7410])
hold on
plot(t_event,x_event(2),'MarkerSize',10,'Marker','pentagram','MarkerFaceColor','m','MarkerEdgeColor','k',"Color","white")
ylim([-0.05 0.4])
legend("","","","Event",'Interpreter','latex');
subplot(2,2,3)
plot(t,x(:,3),'LineWidth',2)
hold on
plot(t_event,x_event(3),'MarkerSize',10,'Marker','pentagram','MarkerFaceColor','m','MarkerEdgeColor','k',"Color","white")
grid on; title('$V_{acc}$(t)','Interpreter','latex')
xlabel('$Time \,\, [s]$','Interpreter','latex'); ylabel('$V_{acc} \,\, [m^3]$','Interpreter','latex');
ylim([8.3e-3 8.9e-3])
legend("","Event",'Interpreter','latex');
subplot(2,2,4)
plot(t,x(:,4),'LineWidth',2)
hold on
plot(t_event,x_event(4),'MarkerSize',10,'Marker','pentagram','MarkerFaceColor','m','MarkerEdgeColor','k',"Color","white")
grid on; title('$V_{tank}$(t)','Interpreter','latex')
xlabel('$Time \,\, [s]$','Interpreter','latex'); ylabel('$V_{tank} \,\, [m^3]$','Interpreter','latex');
ylim([0.9e-3 1.4e-3])
legend("","Event",'Interpreter','latex');
figure()
subplot(1,2,1)
plot(t,rad2deg(param(:,11)),"LineWidth",2)
grid on; title('$\alpha$(t)','Interpreter','latex')
xlabel('$Time \,\, [s]$','Interpreter','latex'); ylabel('$\alpha \,\, [deg]$','Interpreter','latex');
ylim([-20 390]);
subplot(1,2,2)
plot(t,param(:,12),"LineWidth",2)
grid on; title('$A_{dist}$(t)','Interpreter','latex')
xlabel('$Time \,\, [s]$','Interpreter','latex'); ylabel('$A_{dist} \,\, [m^2]$','Interpreter','latex');
ylim([-0.1e-5 2.1e-5]);

figure()
subplot(1,2,1)
plot(t,param(:,1),"LineWidth",2)
hold on
plot(t,param(:,2),"LineWidth",2,"LineStyle","-.")
hold on
plot(t,param(:,3),"LineWidth",2,"LineStyle",":")
hold on
plot(t,param(:,4),"LineWidth",2,"LineStyle","--")
grid on
legend('$P_A$','$P_2$', '$P_3$', '$P_4$','Interpreter','latex')
title('Pressure behaviour in the upper branch','Interpreter','latex')
xlabel('$Time \,\, [s]$','Interpreter','latex'); ylabel('$P \,\, [Pa]$','Interpreter','latex');
subplot(1,2,2)
plot(t,param(:,5),"LineWidth",2)
hold on
plot(t,param(:,6),"LineWidth",2,"LineStyle","-.")
hold on
plot(t,param(:,7),"LineWidth",2,"LineStyle",":")
hold on
plot(t,param(:,8),"LineWidth",2,"LineStyle","--")
grid on
legend('$P_5$', '$P_6$', '$P_7$', '$P_T$','Interpreter','latex')
title('Pressure behaviour in the lower branch','Interpreter','latex')
xlabel('$Time \,\, [s]$','Interpreter','latex'); ylabel('$P \,\, [Pa]$','Interpreter','latex');

figure()
plot(t,param(:,9),"LineWidth",2)
hold on
plot(t,param(:,10),"LineWidth",2)
grid on
legend('$Q_1$', '$Q_2$','Interpreter','latex')
title('Volumetric flow rate','Interpreter','latex')
xlabel('$Time \,\, [s]$','Interpreter','latex'); ylabel('$Q \,\, [m^3/s]$','Interpreter','latex');

figure()
plot(t,param(:,13),"LineWidth",2)
grid on
title('Actuator force','Interpreter','latex')
xlabel('$Time \,\, [s]$','Interpreter','latex'); ylabel('$F \,\, [N]$','Interpreter','latex');
ylim([0,2.6e4])

%% Exercise 3
clearvars; close all; clc;

% Data:
R_1 = 1000;         %[Ohm] Resistor 1
R_2 = 100;          %[Ohm] Resistor 2
L = 0.001;          %[H]   Inductance
C = 0.001;          %[F]   Capacitance
f = 5;              %[Hz]  Voltage source frequency
Vc_0 = 1;           %[V]   Capacitor initial voltage
tspan = [0 10];     %[s]   Time vector

% x = [Vc, Vc_dot]  % State
x0 = [Vc_0; 0];     % Initial conditions

% Case 1: switch closed
odefun_1 = @(t,x) [x(2);
                   (-x(1)-x(2)*(L/R_1 + C*R_2))/(L*C+L*C*R_2/R_1)];
[t_1,x_1] = ode15s(odefun_1,tspan,x0);

% Case 2: voltage source instead of switch
v = @(t) sin(2*pi*f*t)*atan(t);
v_dot = @(t) cos(2*pi*f*t)*(2*pi*f)*atan(t) + sin(2*pi*f*t)/(1+t^2);

odefun_2 = @(t,x) [x(2);...
                   (-x(1)-x(2)*(L/R_1 + C*R_2) + v(t) + v_dot(t)*L/R_1)/(L*C+L*C*R_2/R_1)];

[t_2,x_2] = ode15s(odefun_2,tspan,x0);

% Plots:
figure()
plot(t_1,x_1(:,1),"LineWidth",2)
grid on
title('$V_c$ (with switch closed)','Interpreter','latex')
xlabel('$Time \,\, [s]$','Interpreter','latex'); ylabel('$V_c \,\, [V]$','Interpreter','latex');

figure()
plot(t_2,x_2(:,1),"LineWidth",2)
grid on
title('$V_c$ (with voltage source)','Interpreter','latex')
xlabel('$Time \,\, [s]$','Interpreter','latex'); ylabel('$V_c \,\, [V]$','Interpreter','latex');

%% Exercise 4
clearvars; close all; clc;

% Data:
T_0 = 20;                                     %[°C] Initial temperature (all nodes)
Ti_F = 1000;                                  %[°C] Final temperature (internal node)
To_F = 20;                                    %[°C] Final temperature (outer node)
tspan = [0, 60];                              %[s] Time vector
A = 0.5;                                      %[m^2] Layers area           
l =     [0.001   0.1  0  0.05  0.001];        %[m] Layers thickness vector
k =     [   80   413  0    20     80];        %[W/m/°C] Thermal conductivity
c =     [    0   387  0   470      0];        %[J/kg/°C]  Specific heat
rho =   [    0  8790  0  4500      0];        %[kg/m^3] Layers density

C = A*rho.*c.*l;                              %[J/°C] Thermal capacitance
R = l./k./A; 
R(3) = 0.001;

% Internal node temperature profile:
Ti = @(t) (Ti_F-T_0)*ramp(t,0,1)+T_0;

% Case #1: single node for each layer
R_iA = R(1)/2;
R_AB = R(1)/2 + R(2)/2;
R_BC = R(2)/2 + R(3)/2;
R_CD = R(3)/2 + R(4)/2;
R_DE = R(4)/2 + R(5)/2;
R_Eo = R(5)/2;

dTdt_1 = @(t,x) [( (Ti(t)-x(1))/(R_iA + R_AB) - (x(1)-x(2))/(R_BC + R_CD))/C(2);
                 ( (x(1)-x(2)) /(R_BC + R_CD) - (x(2)-To_F)/(R_DE + R_Eo))/C(4)];

% x = [T_B, T_D];       % state vector
x0_1 = T_0*ones(1,2);   % initial conditions
[t_1,T_1] = ode23s(dTdt_1,tspan,x0_1);

% Nodes array:
Nodes_1 = zeros(length(t_1),7);
Nodes_1(:,1) = Ti(t_1);
Nodes_1(:,3) = T_1(:,1); 
Nodes_1(:,5) = T_1(:,2);
Nodes_1(:,7) = To_F*ones(length(To_F),1);
for i = 1:length(t_1)
    Nodes_1(i,2) = Nodes_1(i,1) - R_iA * (Nodes_1(i,1)-Nodes_1(i,3)) / (R_iA + R_AB);
    Nodes_1(i,4) = Nodes_1(i,3) - R_BC * (Nodes_1(i,3)-Nodes_1(i,5)) / (R_BC + R_CD);
    Nodes_1(i,6) = Nodes_1(i,5) - R_DE * (Nodes_1(i,5)-Nodes_1(i,7)) / (R_DE + R_Eo);
end

% Case #2: two internal nodes for layers B and D
R_iA   = R(1)/2;
R_AB1  = R(1)/2 + R(2)/3;
R_B1B2 = R(2)/3;
R_B2C  = R(2)/3 + R(3)/2;  
R_CD1  = R(3)/2 + R(4)/3;
R_D1D2 = R(4)/3;
R_D2E  = R(4)/3 + R(5)/2;
R_Eo   = R(5)/2;

dTdt_2 = @(t,x) [( (Ti(t)-x(1))/(R_iA + R_AB1)  - (x(1)-x(2))/(R_B1B2)       )/(C(2)/2);
                 ( (x(1)-x(2)) /(R_B1B2)        - (x(2)-x(3))/(R_B2C + R_CD1))/(C(2)/2);
                 ( (x(2)-x(3)) /(R_B2C + R_CD1) - (x(3)-x(4))/(R_D1D2)       )/(C(4)/2);
                 ( (x(3)-x(4)) /(R_D1D2)        - (x(4)-To_F)/(R_D2E + R_Eo) )/(C(4)/2) ];

% x = [T_B1, T_B2, T_D1, T_D2];   % state vector
x0_2 = T_0*ones(1,4);             % initial conditions
[t_2,T_2] = ode23s(dTdt_2,tspan,x0_2);

% Nodes array:
Nodes_2 = zeros(length(t_2),9);
Nodes_2(:,1) = Ti(t_2);
Nodes_2(:,3) = T_2(:,1);
Nodes_2(:,4) = T_2(:,2);
Nodes_2(:,6) = T_2(:,3);
Nodes_2(:,7) = T_2(:,4);
Nodes_2(:,9) = To_F*ones(length(To_F),1);
for i = 1:length(t_2)
    Nodes_2(i,2) = Nodes_2(i,1) - R_iA  * (Nodes_2(i,1)-Nodes_2(i,3)) / (R_iA + R_AB1);
    Nodes_2(i,5) = Nodes_2(i,4) - R_B2C * (Nodes_2(i,4)-Nodes_2(i,6)) / (R_B2C + R_CD1);
    Nodes_2(i,8) = Nodes_2(i,7) - R_D2E * (Nodes_2(i,7)-Nodes_2(i,9)) / (R_D2E + R_Eo);
end

% Plots:
figure()
plot(t_1,Nodes_1(:,1),"LineWidth",2,"LineStyle","-");
hold on
plot(t_1,Nodes_1(:,2),"LineWidth",2,"LineStyle","--");
hold on
plot(t_1,Nodes_1(:,3),"LineWidth",2,"LineStyle","-");
hold on
plot(t_1,Nodes_1(:,4),"LineWidth",2,"LineStyle","--");
hold on
plot(t_1,Nodes_1(:,5),"LineWidth",2,"LineStyle","-");
hold on
plot(t_1,Nodes_1(:,6),"LineWidth",2,"LineStyle","-");
hold on
plot(t_1,Nodes_1(:,7),"LineWidth",2,"LineStyle","--");
grid on
legend("$T_i$","$T_A$","$T_B$","$T_C$","$T_D$","$T_E$","$T_0$",'Interpreter','latex')
title('Node temperature profile','Interpreter','latex')
xlabel('Time [s]','Interpreter','latex'); ylabel('Temperature [$^\circ$C]','Interpreter','latex');
ylim([0 1020])

figure()
plot(t_2,Nodes_2(:,1),"LineWidth",2,"LineStyle","-");
hold on
plot(t_2,Nodes_2(:,2),"LineWidth",2,"LineStyle","--");
hold on
plot(t_2,Nodes_2(:,3),"LineWidth",2,"LineStyle","-.");
hold on
plot(t_2,Nodes_2(:,4),"LineWidth",2,"LineStyle","-");
hold on
plot(t_2,Nodes_2(:,5),"LineWidth",2,"LineStyle","--");
hold on
plot(t_2,Nodes_2(:,6),"LineWidth",2,"LineStyle","-.");
hold on
plot(t_2,Nodes_2(:,7),"LineWidth",2,"LineStyle","-");
hold on
plot(t_2,Nodes_2(:,8),"LineWidth",2,"LineStyle","-");
hold on
plot(t_2,Nodes_2(:,9),"LineWidth",2,"LineStyle","--");
grid on
legend("$T_1$","$T_A$","$T_{B1}$","$T_{B2}$","$T_C$","$T_{D1}$","$T_{D2}$","$T_E$","$T_0$",'Interpreter','latex')
title('Node temperature profile','Interpreter','latex')
xlabel('Time [s]','Interpreter','latex'); ylabel('Temperature [$^\circ$C]','Interpreter','latex');
ylim([0 1020])

%% Functions
function [acc,t,x] = acceleration(x0,Data,k,b)
%
% Function that computes the angular acceleration in time from k and b
%
% INPUTS:
% - x0:      initial value 
% - data:    data structure
% - k:       stiffness coefficent
% - b:       viscous friction coefficent
%
% OUTPUTS:
% - acc:     angular acceleration vector
% - t:       time vector
% - x:       result vector
%
odefun = @(t,x,k,b) [x(3);
                     x(4);
                     k*(x(2)-x(1))/Data.J_1;
                     (-k*(x(2)-x(1)) - sign(x(4))*b*(x(4))^2 + Data.T_0)/Data.J_2];
[t,x] = ode45(@(t,x) odefun(t,x,k,b),Data.tspan,x0);

acc = zeros(length(t),2);
for i = 1:length(t)
    acc(i,1) = k*(x(i,2)-x(i,1))/Data.J_1;
    acc(i,2) = (-k*(x(i,2)-x(i,1)) - sign(x(i,4))*b*(x(i,4))^2 + Data.T_0)/Data.J_2;
end
end

function z = ramp(t,t1,t2)
%
% Function describing a unitary ramp command:
% (z=0 for t<t1, z=ramp for t1 < t < t2, z=1 for t>t2)
%
% INPUTS:
% - t:    time (vector or scalar)
% - t1:   start time of ramp section
% - t2:   end time of ramp section
%
% OUTPUTS:
% - z:    command value
%
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

function [xdot,paramout] = pressurization_BE(t,xx,data)
%
% Function for dynamic simulation of the fluid-mechanical system up until the event
%
% INPUTS:
% - t:           time vector
% - xx:          state vector
% - data:        data structure 
%
% OUTPUTS:
% - xdot:        ODE
% - paramout:    output parameters vector: 
%                paramout = [P_A, P_2, P_3, P_4, P_5, P_6, P_7, P_T, Q1, Q2, alpha, A_dist, F];
%
x = xx(1);
v = xx(2);

% max x constraint:
if x < 0
    x = 0;
end
if x <= 0 && v < 0
    x = 0;
    v = 0;
end
if x > data.piston.x_max
    x = data.piston.x_max;
end
if x >= data.piston.x_max && v > 0
    x = data.piston.x_max;
    v = 0;
end

% Alpha and distributor area behaviour wrt z(t):
alpha = 2*acos(1-2*ramp(t,data.t1,data.t2));
A_dist = data.dist.D^2/8*(alpha-sin(alpha));

% Flow rates
Q1 = data.piston.A1*v;
Q2 =  data.piston.A2*v;

% Set a tolerance to avoid dividing by zero
tol = 1e-10;
if A_dist > tol
    v1_dist = Q1*abs(Q1)/(A_dist)^2;
    v2_dist = Q2*abs(Q2)/(A_dist)^2;
else 
    v1_dist = 0;
    v2_dist = 0;
end

% Pressures computation
P_A = data.gas.p0*(data.gas.V0/(data.accum.V-xx(3)))^data.gas.gamma;           
P_2 = P_A - 0.5*data.fluid.rho*Q1*abs(Q1)*(data.accum.k+data.CV.k)/data.B23.A^2;
P_3 = P_2 - 0.5*data.fluid.rho*Q1*abs(Q1)*data.B23.f*data.B23.L/data.B23.D/data.B23.A^2;
P_T = data.tank.p0;
P_7 = P_T + 0.5*data.fluid.rho*Q2*abs(Q2)*data.tank.k/data.B67.A^2;
P_6 = P_7 + 0.5*data.fluid.rho*Q2*abs(Q2)*data.B67.f*data.B67.L/data.B67.D/data.B67.A^2;
P_5 = P_6 + 0.5*data.fluid.rho*data.dist.k*v2_dist;

if A_dist >= tol
    P_4 = P_3 - 0.5*data.fluid.rho*data.dist.k*v1_dist;
else
    P_4 = (P_T*data.piston.A2 + data.piston.F(0))/data.piston.A1;
end

% xdot:
xdot(1,1) = v;
xdot(2,1) = (P_4*data.piston.A1 - P_5*data.piston.A2 - data.piston.F(x))/data.piston.m;
xdot(3,1) = -Q1;
xdot(4,1) = Q2;
paramout = [P_A, P_2, P_3, P_4, P_5, P_6, P_7, P_T, Q1, Q2, alpha, A_dist, data.piston.F(x)];
end

function [value, isterminal, direction] = press_event(~, xx, data)
%
% Event function for ODE
%
% INPUTS:
% - t:       time (not needed)
% - xx:      state vector
% - data:    data struct 
%
% OUTPUTS:
% - value:       value of event function
% - isterminal:  if 1, integration terminates at value = 0
% - direction:   0     -> detects all events
%                +/- 1 -> detects only events with increasing/decreasing
%                         event function
%
value = data.piston.x_max - xx(1) ;
isterminal = 1 ;
direction = -1 ;
end

