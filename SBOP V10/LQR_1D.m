%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/08/22

%% Set up
set_graphics(); 
close all

%% Problem definition 
% Numerical solver definition 
time_distribution = 'Chebyshev';        % Distribution of time intervals
basis = 'Chebyshev';                    % Polynomial basis to be use
n = 10;                                % Polynomial order in the state vector expansion
m = 60;                                % Number of sampling points
L = 2;                                 % Degree of the dynamics 

OptProblem = Problem().DefineSolver(n, basis, m, time_distribution).AddDynamics(length(n), 1, L); 

% Add boundary conditions
S0 = [0 0].';
SF = [1 0].';
OptProblem = OptProblem.AddBoundaryConditions(S0, SF).AddParameters([]);

% Add functions 
OptProblem = OptProblem.AddFunctions(@(initial, final, beta, t0, tf)BoundaryConditions_LQR(initial, final, beta, t0, tf), @(params, beta, t0, tf, tau, s)ControlFunction_LQR(params, beta, t0, tf, tau, s), ...
                                     @(params, beta, t0, tf, s, u)CostFunction_LQR(params, beta, t0, tf, s, u), @(beta, P)LinConstraints_LQR(beta, P), ...
                                     @(params, beta, t0, tf, tau, s, u)NlinConstraints_LQR(params, beta, t0, tf, tau, s, u), ...
                                     @BoundsFunction_LQR, ...
                                     @(params, initial, final)InitialGuess_LQR(params, initial, final));

%% Optimization
% Simple solution    
tic
[C, cost, u, t0, tf, tau, exitflag, output] = sb_solver(OptProblem);
toc 

% Analytical solution
uopt = -12*tau+6;

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