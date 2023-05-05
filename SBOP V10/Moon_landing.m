%% Project: SBOPT %%
% Date: 05/05/23

%% Moon Landing %% 
% This script provides a main interface to solve the moon landing problem %

%% Set up
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 10;                                % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 2;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 1;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 1;           % Dimension of the control vector

% Boundary conditions
S0 = [100; -20];                % Initial conditions
SF = [0; 0];                    % Final conditions

% Problem parameters 
T = 3;                          % Maximum control acceleration [N]
g = 1.6;                        % Maximum acceleration [m/s^2]

problem_params = [g; T];

% Create the problem
OptProblem = Problems.MinTime_1D(S0, SF, L, StateDimension, ControlDimension, problem_params);

%% Optimization
% Simple solution    
tic
[C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
toc 

u = u-g; 

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
plot(tau, C(1,:))
xlabel('Flight time')
ylabel('$Y$ coordinate')
hold on
grid on; 

figure_orbits = figure;
hold on
plot(tau, C(2,:))
xlabel('Flight time')
ylabel('$Y$ velocity')
hold on
grid on;

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(tau, sqrt(dot(u,u,1)), 'k','LineWidth',1)
plot(tau, u, 'LineWidth', 0.3)
yline(T, '--k')
xlabel('Flight time')
ylabel('$\mathbf{a}$')
legend('$\|u\|$','u', '$T_{max}$')
grid on;