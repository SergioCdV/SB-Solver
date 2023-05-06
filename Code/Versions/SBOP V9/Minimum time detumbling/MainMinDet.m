%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/08/22

%% Set up
set_graphics(); 
close all

%% Problem definition 
% Numerical solver definition 
time_distribution = 'Chebyshev';       % Distribution of time intervals
basis = 'Chebyshev';                   % Polynomial basis to be use
n = [20 20 20];                        % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points

I = [1 0 0; 0 2 0; 0 0 5.5];          % Inertia dyadic
T = 10;                               % Maximum torque
alpha = 3;                            % Maximum axis-wise angular acceleration

OptProblem = Problem().DefineSolver(n, basis, m, time_distribution).AddDynamics(3, 3, 1); 

% Add boundary conditions
S0 = [10 2 -1].';
SF = [0 0 0].';
OptProblem = OptProblem.AddBoundaryConditions(S0, SF).AddParameters([T alpha reshape(I, [], 9)]);

% Add functions 
OptProblem = OptProblem.AddFunctions(@(initial, final, beta, t0, tf)BoundaryConditions_minDet(initial, final, beta, t0, tf), @(params, beta, t0, tf, tau, s)ControlFunction_minDet(params, beta, t0, tf, tau, s), ...
                                     @(params, beta, t0, tf, s, u)CostFunction_minDet(params, beta, t0, tf, s, u), @(beta, P)LinConstraints_minDet(beta, P), ...
                                     @(params, beta, t0, tf, tau, s, u)NlinConstraints_minDet(params, beta, t0, tf, tau, s, u), ...
                                     @BoundsFunction_minDet, ...
                                     @(params, initial, final)InitialGuess_minDet(params, initial, final));

%% Optimization
% Simple solution    
tic
[C, cost, u, t0, tf(1), tau, exitflag, output] = sb_solver(OptProblem);
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