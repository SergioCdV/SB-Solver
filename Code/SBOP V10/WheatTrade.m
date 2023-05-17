%% Project: SBOPT %%
% Date: 01/08/22

%% Planar transfer %% 
% This script provides a main interface to solve in-plane low-thrust transfers in polar coordinates %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = [10; 10];                                % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 1;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 2;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 1;           % Dimension of the control vector

% Boundary conditions
S0 = [1; 1];                % Initial conditions
SF = [0; 0];                    % Final conditions

% Problem parameters 
Tf = 10;              % Planning horizon
Msell = 1;            % Minimum sell price
Mbuy = 20;            % Maximum buy price
R = 0;                % Constant interest rate

problem_params = [Tf; Msell; Mbuy; R];

% Create the problem
OptProblem = Problems.WheatTrading(S0, SF, L, StateDimension, ControlDimension, problem_params);

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
% Radial plot 
hold on
plot(tau, C(1:2,:), 'LineWidth', 0.3)
xlabel('Flight time')
ylabel('$\mathbf{s}$')
grid on;

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(tau, u, 'LineWidth', 0.3)
xlabel('Flight time')
ylabel('$\mathbf{a}$')
legend('$a$','$a_\rho$','$a_\theta$','$a_z$')
grid on;
