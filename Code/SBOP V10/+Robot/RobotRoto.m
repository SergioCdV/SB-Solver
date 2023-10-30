%% Project: SBOPT %%
% Date: 24/10/23

%% 3-DoF optimization of the chaser motion %% 
% This script provides a main interface to solve the general 3-DoF attitude problem for the robot %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 7;                                 % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

Lc = 1;                         % Characteristic length [m]
Tc = 300;                       % Characteristic time [s]
Tmax = 0.004;                   % Maximum available torque [Nm]
omega_max = 0.1;                % Maximum angular velocity [rad/s]

%% Problem definition 
% Add attitude boundary conditions
q0 = [zeros(3,1); 1];                                              % Initial relative quaternion
omega_0 = zeros(4,1);                                              % Initial relative angular velocity of the chaser [rad/s]
omega_0 = 0.5 * QuaternionAlgebra.right_isoclinic(omega_0) * q0;   % Quaternion kinematics
S0 = [q0; omega_0];                                                % Initial conditions

qf = [0.5;0.5;0.5;0.5];                                            % Final relative quaternion (null)
omega_f = zeros(4,1);                                              % Final relative angular velocity [rad/s]
omega_f = 0.5 * QuaternionAlgebra.right_isoclinic(omega_f) * qf;   % Quaternion kinematics
SF = [qf; omega_f];                                                % Final conditions

%% Create the problem
% Linear problem data
params(1) = Tc;                  % TOF [s]
params(2) = Lc;                  % Maximum length [m]
params(3) = Tmax;                % Maximum control authority [Nm]
params(4) = omega_max;           % Maximum rotational speed [rad/s]

% Inertia tensor of the chaser [kg m^2]
params(5:13) = diag([0.001 0.06 0.06]);    

L = 2;                           % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 4;              % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;            % Dimension of the control vector

OptProblem = Problems.RobotRotoBerthing(S0, SF, L, StateDimension, ControlDimension, params);

%% Optimization
% Simple solution    
tic
[C, dV, u, t0, tf, tau, exitflag, output, P0] = solver.solve(OptProblem);
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

% Angular velocity field
omega = 2 * [C(4,:).*C(5,:)+C(3,:).*C(6,:)-C(2,:).*C(7,:)-C(1,:).*C(8,:); ...
            -C(3,:).*C(5,:)+C(4,:).*C(6,:)+C(1,:).*C(7,:)-C(2,:).*C(8,:); ...
             C(2,:).*C(5,:)-C(1,:).*C(6,:)+C(4,:).*C(7,:)-C(3,:).*C(8,:)];

% Dimensions
% C(1:3,:) = C(1:3,:) * Lc;
% tau = tau * 1000;

%% Plots
% State representation
figure
hold on
xlabel('$t$')
ylabel('$\mathbf{q}$')
plot(tau, C(1:4,:));
legend('$q_1$', '$q_2$', '$q_3$')
hold off
legend('off')
grid on;
xlim([0 tau(end)])

figure
hold on
xlabel('$t$')
ylabel('$\mathbf{\omega}$')
plot(tau, omega);
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
yline(Tmax, 'k--')
xlabel('$t$')
ylabel('$\mathbf{\tau}$')
legend('$\tau_x$', '$\tau_y$', '$\tau_z$', '$\|\mathbf{\tau}\|_2$', '$\tau_{max}$');
grid on;
xlim([0 tau(end)])

%% Matlab CSV
csvwrite("test_1.csv", [C; tau])