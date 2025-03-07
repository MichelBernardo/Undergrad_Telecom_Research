clear all
close all
clc

% Minimum example for fmincon

% Number of OFDM subcarriers
N = 4;

% Initial solution required by the solver. It should comply with the
% problem constraints
x0 = ones(N,1);

% Models in matrix form the total power constraint. I am assuming that the
% total available power is 10 W. This constraint models: p1 + ... + pn <= P
A = ones(1,N);
b = 10;

% Define the handle to the the objective funtion to be passed to the solver
fun = @calcCapacity;

% Define the upper/lower bounds
lb = 0.1* ones(N,1);
ub = inf(N,1);

% Calls the solver
[x, fval, exitflag] = fmincon(fun,x0,A,b,[],[],lb,ub);

% Optimal solution
x

% Objective of the optimal solution
fval

% Flag output from the solver giving the status of the solution. Sometimes
% it warns about problems when searching for solution. In general, when it
% assumes the value 1 everything is fine
exitflag
