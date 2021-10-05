%% Modeling and Simulation of Aerospace Systems
%
% Understanding roundoff error
%
% F. Topputo, Politecnico di Milano

%% 1 - Statement of the problem
%
% We have to sum terms having different order of magnitude, e.g.,
%
% x(t*+h) = 1 + 10^(-3) + 10^(-6) + 10^(-9)
% 
% To see how this works, suppose we perform the following sum
% 
% x(t*+h) = 1.0000000000000000 + 
%           0.0011111111111111 + 
%           0.0000022222222222 + 
%           0.0000000033333333 =
%           ------------------
%           1.0011133366666666

%% 2 - Single precision
%
% x(t*+h) = 1.0000000XXXXXXXXX + 
%           0.0011111XXXXXXXXX + 
%           0.0000023XXXXXXXXX + 
%           0.0000000XXXXXXXXX =
%           ------------------
%           1.0011134XXXXXXXXX
% 
% Where 'X' marks the information lost. 
% 
% Note that the error appears at the 7th digit. 

sum_sp = single(1) + single(0.0011111111111111) + single(0.0000022222222222) + single(0.0000000033333333)

%% 2 - Double precision
%
% x(t*+h) = 1.0000000XXXXXXXXX + 
%           0.001111111111111X + 
%           0.000002222222222X + 
%           0.000000003333333X =
%           ------------------
%           1.001113336666667X
% 
% Where 'X' marks the information lost. 
% 
% Note that the error appears at the 15th digit. 

sum_dp = double(1) + double(0.0011111111111111) + double(0.0000022222222222) + double(0.0000000033333333)
