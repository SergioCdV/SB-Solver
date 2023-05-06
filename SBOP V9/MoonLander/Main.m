%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/08/22

%% Set up
set_graphics(); 
close all

%% Problem definition 
% Numerical solver definition 
time_distribution = 'Chebyshev';        % Distribution of time intervals
basis = 'Chebyshev';                    % Polynomial basis to be use
n = 10;                                 % Polynomial order in the state vector expansion
m = 100;                                % Number of sampling points
L = 2;                                  % Degree of the dynamics 

StateDimension = 1; 
ControlDimension = 1; 
OptProblem = Problem().DefineSolver(n, basis, m, time_distribution).AddDynamics(StateDimension, ControlDimension, L); 

% Boundary conditions
S0 = [100; -20];
SF = [0; 0];

% Problem parameters 
T = 3;                   % Maximum acceleration 
g = 1.6;                 % Maximum acceleration 

% Add boundary conditions
OptProblem = OptProblem.AddBoundaryConditions(S0, SF).AddParameters([g; T]);

% Add functions 
OptProblem = OptProblem.AddFunctions(@(initial, final, beta, t0, tf)BoundaryConditionsML(initial, final, beta, t0, tf), @(params, beta, t0, tf, tau, s)ControlFunctionML(params, beta, t0, tf, tau, s), ...
                                     @(params, beta, t0, tf, s, u)CostFunctionML(params, beta, t0, tf, s, u), @(beta, P)LinConstraintsML(beta, P), ...
                                     @(params, beta, t0, tf, tau, s, u)NlinConstraintsML(params, beta, t0, tf, tau, s, u), ...
                                     @BoundsFunctionML, ...
                                     @(params, initial, final)InitialGuessML(params, initial, final));

%% Optimization
% Simple solution    
tic
[C, dV, u, t0, tf, tau, exitflag, output] = sb_solver(OptProblem);
toc 

u = u-g; 

% Average results 
iter = 0; 
time = zeros(1,iter);
setup.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, dV, u, t0, tf, tau, exitflag, output] = sb_solver(OptProblem);
    time(i) = toc;
end

time = mean(time);

%% Plots
% Descend representation
figure_orbits = figure;
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