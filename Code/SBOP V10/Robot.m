%% Project: SBOPT %%
% Date: 01/08/22

%% LQR 1D %% 
% This script provides a main interface to solve the 1D LQR problem %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Bernstein';                    % Polynomial basis to be use
time_distribution = 'Bernstein';        % Distribution of time intervals
n = 10;                                 % Polynomial order in the state vector expansion
m = 300;                                % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

Lc = 0.2;   % Characteristic length [m]
Tc = 30;    % Characteristic time [s]

%% Problem definition 
L = 2;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 1;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 1;           % Dimension of the control vector

% Add boundary conditions
S0 = [0 0].';                   % Initial conditions
SF = [1 0].';                   % Final conditions

% Create the problem
params(1) = Tc;
OptProblem = Problems.Robot(S0, SF, L, StateDimension, ControlDimension, params);

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

% Analytical solution
uopt = -12*tau+6;

C(1,:) = C(1,:) * Lc;
tau = tau * 1000;

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
plot(tau, uopt, 'ok');
xlabel('$t$')
ylabel('$\mathbf{u}$')
grid on;

%% Matlab CSV
csvwrite("test_1.csv", [C; tau])