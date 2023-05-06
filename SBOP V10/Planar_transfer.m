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
n = 8;                                % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 2;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 2;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 2;           % Dimension of the control vector

% Boundary conditions
S0 = [1; 0; 0; 1];              % Initial conditions
SF = [1.1; 0; 0; sqrt(1/1.1)];  % Final conditions

% Problem parameters 
T = 1e-3;              % Maximum acceleration 
mu = 1;                % Gravitational parameter

problem_params = [mu; T];

% Create the problem
OptProblem = Problems.Planar_transfer(S0, SF, L, StateDimension, ControlDimension, problem_params);

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
polarplot(C(2,:), C(1,:)) 

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(tau, T * sqrt(dot(u,u,1)), 'k','LineWidth',1)
plot(tau, T * u, 'LineWidth', 0.3)
yline(T, '--k')
xlabel('Flight time')
ylabel('$\mathbf{a}$')
legend('$a$','$a_\rho$','$a_\theta$','$a_z$')
grid on;

figure 
hold on
plot(tau, rad2deg(atan2(u(2,:),u(1,:)))); 
hold off 
grid on;
xlabel('Time')
ylabel('$\theta$')
title('Thrust in-plane angle')
