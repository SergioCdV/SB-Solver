%% Project: SBOPT %%
% Date: 12/05/23

%% Dryfus %% 
% This script provides a main interface to solve Dryfus' problem %

%% Set up
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                   % Polynomial basis to be use. Alternatively: Legendre, Bernestein, Orthogonal Bernstein
time_distribution = 'Legendre';       % Distribution of time intervals. Alternatively: Bernstein, Orthogonal Bernstein, Chebsyhev, Legendre, Linear, Newton-Cotes, Normal, Random, Trapezoidal
n = [10; 10];                           % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 1;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 2;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 1;           % Dimension of the control vector

% Boundary conditions
S0 = [0;1];                    % Initial conditions
SF = [0;0];                    % Final conditions

% Create the problem
OptProblem = Problems.RossTest2(S0, SF, L, StateDimension, ControlDimension);
    
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
scatter(tau, [2/5-64/5*(2+tau).^(-5); 4./(2+tau).^2])
plot(tau, C(1:2,:))
xlabel('Flight time')
ylabel('$X$ coordinate')
hold off
grid on; 

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(tau, sqrt(dot(u,u,1)), 'k','LineWidth',1)
plot(tau, u, 'LineWidth', 0.3)
scatter(tau, -8./(2+tau).^3)
xlabel('Flight time')
ylabel('$\mathbf{u}$')
grid on;