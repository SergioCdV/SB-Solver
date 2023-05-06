%% Project: SBOPT %%
% Date: 01/08/22

%% Hyper Sensitive %% 
% This script provides a main interface to solve the hyper sensitive problem %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 60;                                % Polynomial order in the state vector expansion
m = 200;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 1;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 1;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 1;           % Dimension of the control vector

% Add boundary conditions
S0 = 1.5;       % Initial conditions
SF = 1;         % Final conditions

T = 50;

% Create the problem
OptProblem = Problems.HyperSensitive(S0, SF, L, StateDimension, ControlDimension, T);

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

%% Plots
% State representation
figure_orbits = figure;
hold on
xlabel('$t$')
ylabel('$\mathbf{s}$')
plot(tau, C(1:2,:));
hold off
legend('off')
grid on;

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(tau, u, 'LineWidth', 0.3)
xlabel('$t$')
ylabel('$\mathbf{u}$')
grid on;