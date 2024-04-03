%% Project: SBOPT %%
% Date: 09/03/24

%% 3D low-thrust transfer in the co-orbital CR3BP problem %% 
% This script provides a main interface to solve 3D low-thrust transfers in the co-orbital CR3BP %

%% Set up 
close all
clear

%% Input data
load COOMOT_2024.mat



%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 10;                                % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
Ld = 2;                         % Degree of the dynamics (maximum derivative order of the ODE system)
ControlDimension = 3;           % Dimension of the control vector

% Initial and final clocks 
t0 = 0;                         % Initial clock 
tf = 2 * pi;                    % Final clock
    
% Initial and final conditions 
rt = target_orbit.Trajectory;
ct = chaser_orbit.Trajectory;

% Initial spiral conditions 
c2 = cn(2);
alpha(1) = (c2-2+sqrt(9*c2^2-8*c2))/2;                  % Characteristic root
alpha(2) = (c2-2-sqrt(9*c2^2-8*c2))/2;                  % Characteristic root
lambda = sqrt(alpha(1));                                % Hyperbolic unstable eigenvalue
omega(1) = sqrt(-alpha(2));                             % Hyperbolic stable eigenvalue
omega(2) = sqrt(c2);                                    % Center manifold eigenvalue
k = (omega(1)^2+1+2*c2)/(2*omega(1));                   % Contraint on the planar amplitude   
    
% Mission parameters 
T = 0.5e-3;                     % Maximum acceleration 
T = T/Gamma;                    % Normalized acceleration

%% Optimization
% Preallocation 
iter = 100; 
time = zeros(2,iter); 
cost = zeros(2,iter);
results = zeros(2,iter); 

% Average results 
for i = 1:iter
    % Initial conditions on the chaser orbit 
    idx = randi([1 size(ct,1)]);
    S0 = ct(idx,1:6).';
    SF = rt(1,1:6).';     

    % Full problem solving 
    StateDimension = 3;                 % Dimension of the configuration vector. Note the difference with the state vector
    problem_params = [mu; t0; tf; T];

    OptProblem = Problems.CR3BP_transfer(S0, SF, Ld, StateDimension, ControlDimension, problem_params);
    tic 
    [C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
    time(1,i) = toc;
    cost(1,i) = dV;
    results(1,i) = 1 * (exitflag ~= -2); 

    % Amplitude problem solving
    StateDimension = 2;                       % Dimension of the configuration vector. Note the difference with the state vector
    S0 = ct(idx,1:6)-rt(1,1:6);               % Relative initial particle

    psi = pi/2;                               % Out-of-plane angle
    phi = atan2(S0(2), -k * S0(1));           % In-plane angle
    
    A = [ abs(S0(1) / cos(phi)); S0(3) / sin(psi)];
    S0 = [A; -(S0(4)-A(1)*omega(1)*sin(phi)) / cos(phi); (S0(6)-A(2)*omega(2)*cos(psi)) / sin(psi)];
    SF = zeros(Ld*StateDimension,1);                         % Final particle

    problem_params = [mu; t0; tf; omega.'; k; psi; phi; T];

    OptProblem = Problems.IntegralsTransferCR3BP(S0, SF, Ld, StateDimension, ControlDimension, problem_params);

    tic 
    [C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
    time(2,i) = toc;
    cost(2,i) = dV;
    results(2,i) = 1 * (exitflag ~= -2); 
end


cost = cost * Vc;

%% Plots
figure 
plot(1:i, results(:,1:i))
legend('CR3BP', '$H_2$')
xlabel('Test case')
ylabel('Convergence rate')
grid on;

figure 
plot(1:i, time(:,1:i))
legend('CR3BP', '$H_2$')
xlabel('Test case')
ylabel('Computational cost [s]')
grid on;

figure 
plot(1:i, cost(:,1:i) / Vc)
legend('CR3BP', '$H_2$')
xlabel('Test case')
ylabel('$\Delta V$')
grid on;

%% Save results in a common file 
save COOMOT_2024_comparison_results.mat