%% Project: SBOPT %%
% Date: 12/05/23

%% Breakwell %% 
% This script provides a main interface to solve Breakwell' problem %

%% Set up
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                   % Polynomial basis to be use. Alternatively: Legendre, Bernestein, Orthogonal Bernstein
time_distribution = 'Legendre';       % Distribution of time intervals. Alternatively: Bernstein, Orthogonal Bernstein, Chebsyhev, Legendre, Linear, Newton-Cotes, Normal, Random, Trapezoidal
n = 10;                               % Polynomial order in the state vector expansion
m = 100;                              % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 2;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 1;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 1;           % Dimension of the control vector

% Boundary conditions
S0 = [0; 1];                    % Initial conditions
SF = [0;-1];                    % Final conditions

% Parameters 
TOF = 1;                       % Final time of flight
l = 0.1;                       % State constraint

params = [TOF; l];

% Create the problem
OptProblem = Problems.Breakwell(S0, SF, L, StateDimension, ControlDimension, params);
    
%% Optimization
% Simple solution    
tic
[C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
toc 

% Average results 
iter = 25; 
time = zeros(1,iter);
setup.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
    time(i) = toc;
end

time = mean(time);

% Analytical solution 
t1 = 3*l;
t2 = 1 - 3*l;
uopt = zeros(1,length(tau)); 
uopt(tau < t1) = -2/(3*l)*(1-tau(tau < t1)/(3*l));
uopt(tau > t2) = -2/(3*l)*(1-(1-tau(tau > t2))/(3*l));

%% Plots
% State representation
figure;
hold on
plot(tau, C(1:2,:))
yline(l, 'k--')
xlabel('$t$')
ylabel('$\mathbf{s}$ ')
legend('$s$', '$\dot{s}$', '$s = l$')
hold off
grid on; 
xlim([0 tf])
ylim([min(C(2,:)) max(C(2,:))])

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(tau, uopt, 'ok');
plot(tau, u)

xlabel('$t$')
ylabel('$\mathbf{u}$')
legend('$u^*$','SBOPT')
grid on;