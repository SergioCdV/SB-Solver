%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/08/22

%% Set up
set_graphics(); 
close all

%% Problem definition 
% Numerical solver definition 
time_distribution = 'Chebyshev';        % Distribution of time intervals
basis = 'Chebyshev';                    % Polynomial basis to be use
n = 60;                                % Polynomial order in the state vector expansion
m = 200;                                % Number of sampling points
L = 1;                                 % Degree of the dynamics 

OptProblem = Problem().DefineSolver(n, basis, m, time_distribution).AddDynamics(length(n), 1, L); 

% Add boundary conditions
S0 = 1.5;
SF = 1;
T = 50;
OptProblem = OptProblem.AddBoundaryConditions(S0, SF).AddParameters([T]);

% Add functions 
OptProblem = OptProblem.AddFunctions(@(initial, final, beta, t0, tf)BoundaryConditions_HS(initial, final, beta, t0, tf), @(params, beta, t0, tf, tau, s)ControlFunction_HS(params, beta, t0, tf, tau, s), ...
                                     @(params, beta, t0, tf, s, u)CostFunction_HS(params, beta, t0, tf, s, u), @(beta, P)LinConstraints_HS(beta, P), ...
                                     @(params, beta, t0, tf, tau, s, u)NlinConstraints_HS(params, beta, t0, tf, tau, s, u), ...
                                     @BoundsFunction_HS, ...
                                     @(params, initial, final)InitialGuess_HS(params, initial, final));

%% Optimization
% Simple solution    
tic
[C, cost, u, t0, tf, tau, exitflag, output] = sb_solver(OptProblem);
toc 

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