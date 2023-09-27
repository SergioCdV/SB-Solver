%% Project: SBOPT %%
% Date: 01/08/22

%% Poitning with torquers %% 
% This script provides a main interface to solve the magnetic pointing problem %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                   % Polynomial basis to be use
time_distribution = 'Legendre';       % Distribution of time intervals
n = [15 15 15];                        % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 1;                                % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 3;                   % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;                 % Dimension of the control vector

% Problem parameters
mu = 3.986e14;                        % Gravitational parameter of the Earth
r = 6900e3;                           % Orbital radius 
omega = sqrt(mu/r^3);                 % Angular velocity of the anomaly
i = deg2rad( 98 );                    % Orbital inclination
theta0 = 0;                           % Reference anomaly

Bref = 7.83e15 / r^3;                          % Reference dipole moment
I = [1e-2 0 0; 0 2e-2 0; 0 0 5.5e-2];          % Inertia dyadic
T = [0.2; 0.2; 0.2];                           % Maximum moment per axis
alpha = deg2rad( 90 );                        % Maximum axis-wise angular velocity

params = [omega i theta0 Bref reshape(I, 1, []) T.' alpha].';

% Add boundary conditions
S0 = [deg2rad( 80 ) deg2rad( -5 ) deg2rad( 3 )].';
SF = zeros(3,1);

% Create the problem
OptProblem = Problems.MagneticDetumbling(S0, SF, L, StateDimension, ControlDimension, params);

%% Optimization
% Simple solution    
tic
[C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
toc 

% Magnetic field
B = params(4) * [sin( params(2) ) * cos(tau); repmat( -cos( params(2) ), 1, length(tau) );  2 * sin( params(2) ) * cos(tau)];

% Time law
tau = tau / params(1);

% Average results 
iter = 0; 
time = zeros(1,iter);
setup.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
    time(i) = toc;
end

%% Plots
% State representation
figure_orbits = figure;
hold on
xlabel('$t$')
ylabel('$\mathbf{s}$')
plot(tau, C(1:3,:));
yline(alpha, 'k--');
hold off
legend('off')
grid on;

% Control plot
% Cross-product law 
for i = 1:length(tau)
    u(:,i) = cross(B(:,i), u(:,i)) / norm(B(:,i))^2;
end

figure_propulsion = figure;
hold on
plot(tau, u, 'LineWidth', 0.3)
yline(T, 'k--');
xlabel('$t$')
ylabel('$\mathbf{u}$')
grid on;