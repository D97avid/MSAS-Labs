% Modeling and Simulation of Aerospace Systems (2021/2022)
% Assignment # 1
% Author: David Reina

%% Exercise 1
clearvars; close all; clc;
f = @(x) cos(x) - x;

%% Exercise 2
clearvars; close all; clc;
f = @(x) [x(1)^2 - x(1) - x(2), x(1)^2/16 + x(2)^2 - 1]';

%% Exercise 3
clearvars; close all; clc;
x0 = 0.5;
x_analytic = @(t) t^2 + 2*t + 1 - 0.5*e^t;

%% Exercise 4
clearvars; close all; clc;
A = @(aa) [0 1; -1 2*cos(aa)];
alpha = linspace(0,2*pi,1e3);
lambda = zeros(size(alpha));

for i = 1:length(alpha)
    lambda(i,:) = eig(A(alpha(i)))';
end


%% Exercise 5
clearvars; close all; clc;


%% Exercise 6
clearvars; close all; clc;


%% Exercise 7
clearvars; close all; clc;
% NB B is a STIFF MATRIX

B = [-180.5, 219.5; 179.5, -220.5];
x0 = [1,1]';


%% Functions

