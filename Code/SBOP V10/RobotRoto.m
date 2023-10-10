%% Project: SBOPT %%
% Date: 01/08/22

%% 6-DoF optimization of the chaser motion %% 
% This script provides a main interface to solve the general 6-DoF problem for the robot %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Bernstein';                    % Polynomial basis to be use
time_distribution = 'Bernstein';        % Distribution of time intervals
n = 7;                                 % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

Lc = 1;                         % Characteristic length [m]
Tc = 300;                  % Characteristic time [s]
Tmax = 0.004;                      % Maximum available torque [Nm]

%% Problem definition 
% Target orbital elements
Orbit_t = [7011e3 0.05 deg2rad(190) deg2rad(98) 0 0];

mu = 3.986e14;                  % Gravitational parameter of the Earth
n = sqrt(mu/Orbit_t(1)^3);      % Mean motion
K = floor(Tc/(2*pi/n));         % Number of complete revolutions
dt = Tc - K * (2*pi/n);         % Elapsed time in the last revolution [s]

elements = Orbit_t; 
nu_0 = OrbitalDynamics.kepler(elements);                        % Initial true anomaly [rad]
elements(6) = elements(6) + n * dt;                             % Final mean anomaly
nu_f = 2*pi*K + OrbitalDynamics.kepler(elements);               % Final true anomaly [rad]
h = sqrt(mu * Orbit_t(1) * (1-Orbit_t(2)^2));                   % Target angular momentum

nu = linspace(nu_0, nu_f, m+1);                                   % Anomaly space

% Add attitude boundary conditions
sigma_0 = [zeros(3,1); 1];                                              % Initial relative MRP
omega_0 = zeros(4,1);                                                   % Initial relative angular velocity of the chaser [rad/s]
omega_0 = 0.5 * QuaternionAlgebra.right_isoclinic(sigma_0) * omega_0;   % Quaternion kinematics
S0 = [sigma_0; omega_0];                                                % Initial conditions

sigma_f = [0.5;0.5;0.5;0.5];                                                    % Final relative quaternion (null)
omega_f = zeros(4,1);                                                   % Final relative angular velocity [rad/s]
omega_0 = 0.5 * QuaternionAlgebra.right_isoclinic(sigma_f) * omega_f;   % Quaternion kinematics
SF = [sigma_f; omega_f];                                                % Final conditions

%% Create the problem
% Linear problem data
params(1) = Tc;                  % TOF 
params(2) = Lc;                  % Maximum length
params(3) = Tmax;                % Maximum control authority 

params(4) = mu;                  % Gauss constant
params(5) = Orbit_t(2);          % Target orbital eccentricity
params(6) = h;                   % Angular momentum magnitude
params(7) = nu_0;                % Initial true anomaly [rad]
params(8) = nu_f;                % Final true anomaly [rad]

params(9:17) = diag([0.001 0.06 0.06]);    % Inertia tensor of the chaser [kg m^2]

params = [params nu];            % Final parameter vector

L = 2;                           % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 4;              % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;            % Dimension of the control vector

OptProblem = Problems.RobotRotoBerthing(S0, SF, L, StateDimension, ControlDimension, params);

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