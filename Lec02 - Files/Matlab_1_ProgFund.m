%% Modeling and Simulation of Aerospace Systems
% Some notes on Matlab: Review of programming fundamentals
%
% This is a lecture to recall some basic notions of Matlab programming
%
% F. Topputo, Politecnico di Milano

%% 1 - Obtaining and insatalling Matlab with Polimi license (see slides)
%
% Go to:
%%
% # polimi.it
% # Online services
% # University ICT services
% # Software download
% # Students
% # Matlab
%
% ... and follow the steps therein

%% 2 - Matlab environment (see slides)

%%
% 
% * Command window
% * Editor
% * Workspace
% * Current folder
% * Command history
% * The help
% * ...

%% 3 - Matlab classes

%% 3.1 - Floating-point numbers
clear all
a = 4.52           ; % double precision floating-point number (default)
b = single(4.52)   ; % single precision floating-point number
Af = rand(10,10)   ; % full matrix
As = sparse(10,10) ; % sparse matrix, see the attribute
whos

%% 3.2 - Integers
clear all
c = int8(2.3)   ; % signed, 8-bit integer
d = int8(1000)  ; % check value assigned to d!
maxint8 = intmax('int8')    % largest positive 8-bit integer value
e = int16(1000) ; % signed, 8-bit integer
maxint16 = intmax('int16')   % largest positive 16-bit integer value
f = int32(-10)  ; % signed, 32-bit integer
g = uint32(-10)  ; % unsigned, 32-bit integer, check the value assigned
whos

%% 3.3 - Characters
clear all
h = 'ciao'            ; % simple string
i = char('pippo')     ; % another string of characters
l = strcat([h,' ',i]) ; % horizontal concatenation of strings
whos 

%% 3.4 - Logical
clear all
m = logical(123) ;
n = false        ;
n == m   
whos

%% 3.5 - Function handle
clear all
fun1 = @(x) sin(3*x) + x^2 ; % handle to anonymous function
fun1(3.1)                    % evaluation
fun2 = @sin                ; % handle to a named function
fun2(pi/3)                   % evaluation
whos

%% 3.6 - Tables
clear all
LastName = {'Smith';'Johnson';'Williams';'Jones';'Brown'};
Age = [38;43;38;40;49];
Height = [71;69;64;67;64];
Weight = [176;163;131;133;119];
BloodPressure = [124 93; 109 77; 125 83; 117 75; 122 80];
myTable = table(Age,Height,Weight,BloodPressure,'RowNames',LastName)
whos

%% 3.7 - Structures 
clear all
myVar.mass = 23 ;
myVar.length = 12 ;
myVar.bounds = [12 34] ;
myVar.name = 'pippo' ;
myVar.matrix = rand(2,2) ;
myVar
whos

%% 3.8 - Cell arrays
clear all
myCell = {1, 2, 3;
          'text', rand(5,10,2), {11; 22; 33}} ;
myCell{1,1}, myCell{1,2}, myCell{1,3}
myCell{2,1}, myCell{2,2}, myCell{2,3}
whos

%% 4 - Graphics
x = linspace(0,2*pi,50) ;
fp = plot(x, sin(x), 'b-.', 'Linewidth', 2, 'Marker', 'o') ;
set(gca, 'FontSize', 20) ;
title('Sample figure') ;
xlabel('$x$', 'Interpreter', 'Latex')
ylabel('$\sin x$', 'Interpreter', 'Latex') ;
legend('sin x') ;

%% 5 - Programming in Matlab
%
% 
% * *Scripts*
% * Scripts > Sections
% * *Functions*
% * Functions > Local
% * Functions > Nested
% * Functions > Anonymous

%% 6 - Static vs dynamic allocation
% Exampple of inefficient dynamic allocation
clear all
tic
for i=1:1000
       for j=1:5000
           matrix1(i,j) = i*(j+1);
       end
end
time1 = toc ;
% Static allocation
tic
matrix2 = zeros (1000, 1000) ;
for i=1:1000
       for j=1:5000
           matrix2(i,j) = i*(j+1);
       end
end
time2 = toc ;
fprintf('Dynamic allocation, t = %d s\n', time1) ;
fprintf('Static allocation, t = %d s\n', time2) ;
fprintf('Speedup = %d\n', time1/time2) ;
% Example of vectorization
tic ;
row = 1:5000;
col = 1:1000;
matrix3 = row'*(col+1);
time3 = toc ;
fprintf('Vectorization, t = %d s\n', time3) ;
fprintf('Speedup = %d\n', time2/time3) ;