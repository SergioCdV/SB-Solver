%% Project: SBOPT %%
% Date: 01/08/22

%% Detumbling with torquers %% 
% This script provides a main interface to solve the angular momentum detumbling problem %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = [20 20 20];                        % Polynomial order in the state vector expansion
m = 200;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 1;                                % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 3;                   % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;                 % Dimension of the control vector

% Problem parameters
I = [1 0 0; 0 2 0; 0 0 5.5];          % Inertia dyadic
T = 10;                               % Maximum torque
alpha = 3;                            % Maximum axis-wise angular acceleration

problem_params = [T alpha reshape(I, [], 9)];

% Add boundary conditions
S0 = [10 2 -1].';
SF = [0 0 0].';

% Create the problem
OptProblem = Problems.HyperSensitive(S0, SF, L, StateDimension, ControlDimension, problem_params);

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