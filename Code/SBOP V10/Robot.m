%% Project: SBOPT %%
% Date: 01/08/22

%% LQR 1D %% 
% This script provides a main interface to solve the 1D LQR problem %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 10;                                 % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

Lc = 0.2;                       % Characteristic length [m]
Tc = 3600;                       % Characteristic time [s]
Fmax = 2;                       % Maximum available force [m/s^2]
Tmax = 1;                       % Maximum available torque [Nm]

%% Problem definition 
L = 2;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 7;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 6;           % Dimension of the control vector

% Add boundary conditions
sigma_0 = [zeros(3,1); 1];                                      % Initial relative MRP
omega_0 = zeros(4,1);                                           % Initial relative angular velocity of the chaser [rad/s]
sigma_f = 0.5 * ones(4,1);                                      % Final relative MRP
omega_f = zeros(4,1);                                           % Final relative angular velocity [rad/s]
S0 = [zeros(1,3).'; sigma_0; zeros(1,3).'; omega_0];            % Initial conditions
SF = [zeros(1,3).'; sigma_f; zeros(1,3).'; omega_f];            % Final conditions

% Create the problem
params(1) = Tc;                  % TOF 
params(2) = Fmax;                % Maximum control authority 
params(3) = Lc;                  % Maximum length
params(4:6) = [1 1 1].';         % Final target's docking port position 
params(7:9) = [1 1 1].';         % Final target's docking port angular velocity
params(10:12) = [5 5 10].';      % Final target's docking port center of mass
params(13:21) = diag([1 2 3]);   % Inertia tensor of the chaser [kg m^2]
params(22) = Tmax;               % Maximum torques

OptProblem = Problems.Robot(S0, SF, L, StateDimension, ControlDimension, params);

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

% Dimensions
% C(1:3,:) = C(1:3,:) * Lc;
% tau = tau * 1000;

%% Plots
% State representation
figure
hold on
xlabel('$t$')
ylabel('$\mathbf{r}$')
plot(tau, C(1:3,:));
legend('$x$', '$y$', '$z$')
hold off
legend('off')
grid on;
xlim([0 tau(end)])

figure
hold on
xlabel('$t$')
ylabel('$\mathbf{v}$')
plot(tau, C(8:10,:));
legend('$v_x$', '$v_y$', '$v_z$')
hold off
legend('off')
grid on;
xlim([0 tau(end)])

figure
hold on
xlabel('$t$')
ylabel('$\mathbf{q}$')
plot(tau, C(4:7,:));
legend('$q_1$', '$q_2$', '$q_3$')
hold off
legend('off')
grid on;
xlim([0 tau(end)])

figure
hold on
xlabel('$t$')
ylabel('$\mathbf{\omega}$')
plot(tau, C(11:13,:));
legend('$\omega_x$', '$\omega_y$', '$\omega_z$')
hold off
legend('off')
grid on;
xlim([0 tau(end)])

% Propulsive acceleration plot
figure;
hold on
plot(tau, u(1:3,:), 'LineWidth', 0.3)
plot(tau, sqrt(dot(u(1:3,:),u(1:3,:),1)), 'k');
yline(Fmax, 'k--')
xlabel('$t$')
ylabel('$\mathbf{a}$')
legend('$a_x$', '$a_y$', '$a_z$', '$\|\mathbf{a}\|_2$', '$a_{max}$');
grid on;
xlim([0 tau(end)])

figure;
hold on
plot(tau, u(1:3,:), 'LineWidth', 0.3)
plot(tau, sqrt(dot(u(4:6,:),u(4:6,:),1)), 'k');
yline(Tmax, 'k--')
xlabel('$t$')
ylabel('$\mathbf{\tau}$')
legend('$\tau_x$', '$\tau_y$', '$\tau_z$', '$\|\mathbf{\tau}\|_2$', '$\tau_{max}$');
grid on;
xlim([0 tau(end)])

figure
view(3)
hold on
xlabel('$x$')
ylabel('$y$')
ylabel('$z$')
plot3(C(1,:), C(2,:), C(3,:));
hold off
legend('off')
grid on;

%% Matlab CSV
csvwrite("test_1.csv", [C; tau])