%% Project: SBOPT %%
% Date: 05/05/23

%% Zermelos %% 
% This script provides a main interface to solve Zermelo's problem %

%% Set up
close all
clear

%% Numerical solver definition 
basis = 'Bernstein';                   % Polynomial basis to be use. Alternatively: Legendre, Bernestein, Orthogonal Bernstein
time_distribution = 'Linear';       % Distribution of time intervals. Alternatively: Bernstein, Orthogonal Bernstein, Chebsyhev, Legendre, Linear, Newton-Cotes, Normal, Random, Trapezoidal
n = 30;                                % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 2;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 2;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 2;           % Dimension of the control vector

% Boundary conditions
S0 = [1; 0; 0; 0];              % Initial conditions
SF = [0; 5; 45; 0];             % Final conditions

% Problem parameters 
a = 100;                        % Maximum thrust [N]

% Create the problem
OptProblem = Problems.Lintan_steering(S0, SF, L, StateDimension, ControlDimension, a);
    
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
plot(tau, C(1,:))
xlabel('Flight time')
ylabel('$X$ coordinate')
hold off
grid on; 

figure
hold on
plot(tau, C(2,:))
xlabel('Flight time')
ylabel('$Y$ coordinate')
hold off
grid on;

figure_orbits = figure;
hold on
plot(C(1,:), C(2,:))
xlabel('$X$ coordinate')
ylabel('$Y$ coordinate')
hold on
grid on;

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(tau, u, 'LineWidth', 0.3)
yline(a, '--')
xlabel('Flight time')
ylabel('$\mathbf{u}$')
grid on;

figure_propulsion = figure;
hold on
plot(tau, atan2(u(2,:),u(1,:)), 'LineWidth', 0.3)
xlabel('Flight time')
ylabel('$\alpha$')
grid on;