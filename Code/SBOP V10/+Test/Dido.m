%% Project: SBOPT %%
% Date: 12/05/23

%% Dryfus %% 
% This script provides a main interface to solve Dryfus' problem %

%% Set up
close all
clear

%% Numerical solver definition 
basis = 'Chebyshev';                   % Polynomial basis to be use. Alternatively: Legendre, Bernestein, Orthogonal Bernstein
time_distribution = 'Chebyshev';       % Distribution of time intervals. Alternatively: Bernstein, Orthogonal Bernstein, Chebsyhev, Legendre, Linear, Newton-Cotes, Normal, Random, Trapezoidal
n = [10; 10; 10];                      % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 1;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 3;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;           % Dimension of the control vector

% Boundary conditions
S0 = [0; 0; 0];                 % Initial conditions
SF = [1; 0; 1];                 % Final conditions

% Problem parameters
len = 20;

% Create the problem
OptProblem = Problems.Dido(S0, SF, L, StateDimension, ControlDimension, len);
    
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
xlabel('Flight time')
ylabel('$X$ coordinate')
hold off
grid on; 

figure;
view(3)
hold on
plot3(C(1,:), C(2,:), C(3,:))
xlabel('$X$ coordinate')
ylabel('$Y$ coordinate')
ylabel('$Z$ coordinate')
hold off
grid on; 


% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(tau, sqrt(dot(u,u,1)), 'k','LineWidth',1)
plot(tau, u, 'LineWidth', 0.3)
xlabel('Flight time')
ylabel('$\mathbf{u}$')
grid on;