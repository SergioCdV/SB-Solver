%% Project: SBOPT %%
% Date: 10/11/23

%% Planar transfer, Herman, Conway 1996 %% 
% This script provides a main interface to solve in-plane low-thrust transfers in polar coordinates %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = [7; 7];                          % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 2;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 2;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 2;           % Dimension of the control vector

% Boundary conditions
S0 = [1.1; 0; 0; (1 / sqrt(1.1)) / 1.1];      % Initial conditions
SF = 4 * S0;                                    % Final conditions

% Problem parameters 
TOF = 50;                       % Maximum acceleration 
mu = 1;                         % Gravitational parameter
A = 0.01;                       % Maximum control

params = [mu; TOF; A];

% Create the problem
OptProblem = Problems.HermanConway(S0, SF, L, StateDimension, ControlDimension, params);

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
figure 
hold on
theta = linspace(0, 2*pi, 100); 
x_t = C(1,end) * cos(theta);
y_t = C(1,end) * sin(theta);
plot(C(1,:) .* cos(C(2,:)), C(1,:) .* sin(C(2,:)) ) 
plot(x_t, y_t, 'k--') 
xlabel('$x$')
ylabel('$y$')
grid on;

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(tau, atan2(u(1,:), u(2,:)), 'k');
xlabel('$t$')
ylabel('$\beta$')
grid on;
xlim([0 tf])
