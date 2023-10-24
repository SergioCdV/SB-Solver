%% Project: SBOPT %%
% Date: 01/08/22

%% LQRC 1D %% 
% This script provides a main interface to solve a constrained 1D LQR problem %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Bernstein';                    % Polynomial basis to be use
time_distribution = 'Bernstein';        % Distribution of time intervals
n = 80;                                % Polynomial order in the state vector expansion
m = 300;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 2;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 1;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 1;           % Dimension of the control vector

% Add boundary conditions
S0 = [0 1].';                   % Initial conditions
SF = [0 -1].';                  % Final conditions

% Problem parameters 
A = 1/6;

% Create the problem
OptProblem = Problems.Breakwell(S0, SF, L, StateDimension, ControlDimension, A);

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