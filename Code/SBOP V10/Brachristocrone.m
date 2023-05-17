%% Project: SBOPT %%
% Date: 05/05/23

%% Zermelos %% 
% This script provides a main interface to solve Zermelo's problem %

%% Set up
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                   % Polynomial basis to be use. Alternatively: Legendre, Bernestein, Orthogonal Bernstein
time_distribution = 'Legendre';       % Distribution of time intervals. Alternatively: Bernstein, Orthogonal Bernstein, Chebsyhev, Legendre, Linear, Newton-Cotes, Normal, Random, Trapezoidal
n = [9 9 9];                           % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 1;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 3;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 1;           % Dimension of the control vector

% Boundary conditions
S0 = [0; 0; 0];                 % Initial conditions
SF = [10; 10; 0];                 % Final conditions

% Problem parameters 
g = 9.81;                       % Gravity acceleration

problem_params = g;

% Create the problem
OptProblem = Problems.Brachristocrone(S0, SF, L, StateDimension, ControlDimension, problem_params);
    
%% Optimization
% Simple solution    
tic
[C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
toc 

% Average results 
iter = 0; 
time = zeros(1,iter);
setup.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
    time(i) = toc;
end

time = mean(time);

%% Plots
% Descend representation
figure;
hold on
plot(tau, C(1:3,:))
xlabel('$X$ coordinate')
ylabel('$Y$ coordinate')
hold off
grid on;

figure;
hold on
plot(C(1,:), C(2,:))
xlabel('$X$ coordinate')
ylabel('$Y$ coordinate')
hold off
grid on; 

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(tau, atan2(C(4,:),C(5,:)), 'LineWidth', 0.3)
xlabel('Flight time')
ylabel('$\mathbf{a}$')
grid on;