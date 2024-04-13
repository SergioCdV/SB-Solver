%% Project: SBOPT %%
% Date: 13/04/24

%% Attitude slew with RW %% 
% This script provides a main interface to solve an attitude slew problem usign Reaction Wheels %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 7;                        % Polynomial order in the state vector expansion
m = 50;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 2;                                % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 6;                   % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;                 % Dimension of the control vector

% Problem parameters
t0 = 0;                                % Initial clock 
tf = 120;                              % Final clock
I = [0.6 0.1 0; 0.1 1 0.2; 0 0.2 1] * 1E-2;    % Inertia dyadic
T = 0.0004;                             % Maximum torque [Nm]
omega_max = 0.1;                       % Maximum angular velocity [rad/s]
h_max = 0.002;                         % Maximum angular momentum [Nms]

problem_params = [t0 tf T omega_max h_max reshape(I, [], 9)];

% Add boundary conditions
S0 = [zeros(1,12)].';
SF = [1/3 1/3 1/3 zeros(1,9)].';

% Create the problem
OptProblem = Problems.AttitudeSlew(S0, SF, L, StateDimension, ControlDimension, problem_params);

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

%% Plots
% State representation
figure_orbits = figure;
hold on
xlabel('$t$')
ylabel('$\mathbf{\sigma}$')
plot(tau, C(1:3,:));
hold off
legend('$\sigma_x$', '$\sigma_y$', '$\sigma_z$')
grid on;
yticklabels(strrep(yticklabels, '-', '$-$'));

figure;
hold on
xlabel('$t$')
ylabel('$\mathbf{h}$')
plot(tau, C(4:6,:));
plot(tau, max(abs(C(4:6,:))))
yline(h_max, 'k--');
hold off
legend('$h_x$', '$h_y$', '$h_z$', '$\|\mathbf{h}\|_{\infty}$', '$h_{max}$')
grid on;
yticklabels(strrep(yticklabels, '-', '$-$'));

% Torque plot
figure_propulsion = figure;
hold on
plot(tau, u)
plot(tau, sqrt(dot(u,u,1)))
%yline(T, 'k--');
xlabel('$t$')
ylabel('$\mathbf{\tau}$')
legend('$\tau_x$', '$\tau_y$', '$\tau_z$', '$\|\mathbf{\tau}\|_2$', '$\tau_{max}$')
grid on;
yticklabels(strrep(yticklabels, '-', '$-$'));