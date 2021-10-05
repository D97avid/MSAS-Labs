clear; clc; close all;

%% Stability region
alpha = linspace(0,2*pi,1e3)';
h = zeros(length(alpha),1);
eig_A = [h,h];
lambda = [h,h];

% Integration schemes
A = @(x) [0 1; -1 2*cos(x)];
FE = @(hh,x) eye(2) + hh*A(x);
BE = @(hh,x) (eye(2) - hh*A(x))'; %funziona con la trasposizione, NON con l'inversione inv()

RK4 = @(hh,x) eye(2) + hh/6 * (A(x) + 2*(A(x) + hh/2 * A(x)^2) + 2*(A(x) + hh/2 * A(x)^2 + hh^2/4 * A(x)^3) + (A(x) + hh * A(x)^2 + hh^2/2 * A(x)^3 + hh^3/4 * A(x)^4));

IEX4 = @(hh,x) -1/6 * (eye(2) + hh*A(x)) + 4*(eye(2) + hh/2*A(x))^2 - 27/2*(eye(2) + hh/3*A(x))^3 + 32/3*(eye(2) + hh/4*A(x))^4;


scheme = 'RK4';

switch(scheme)
    case 'FE'
        fun = FE;
    case 'BE'
        fun = BE;
    case 'RK4'
        fun = RK4;
    case 'IEX4'
        fun = IEX4;
end

options = optimoptions('fsolve','Display','none');
%options = optimset('Display','off');
for i = 1:length(alpha)
    h(i) = fsolve(@(hh) (max(abs(eig(fun(hh,alpha(i))))) - 1),5,options);
    %h(i) = fzero(@(hh) (max(abs(eig(fun(hh,alpha(i))))) - 1),5,options);
    lambda(i,:) = (eig(A(alpha(i))))';
end

figure()
plot(real(h.*lambda),imag(h.*lambda))
title("Stability region " + "(" + scheme + ")")
grid on
axis equal
