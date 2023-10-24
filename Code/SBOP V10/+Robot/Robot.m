%% Project: SBOPT %%
% Date: 01/08/22

%% 6-DoF optimization of the chaser motion %% 
% This script provides a main interface to solve the general 6-DoF problem for the robot %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 12;                                % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

Lc = 1;                         % Characteristic length [m]
Tc = 100;                       % Characteristic time [s]
Fmax = 0.5e-1;                  % Maximum available acceleration [m/s^2]

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

% Add linear boundary conditions
r_0 = 10 * rand(1,3);                                           % Initial dimensional position vector [m]
v_0 = 0.1 * rand(1,3);                                          % Initial dimensional velocity vector [m/s]
S0 = [r_0.'; v_0.'];                                            % Initial conditions

%% Final boundary conditions
% TH space transformation 
omega = mu^2 / h^3;                                             % True anomaly angular velocity
k = 1 + Orbit_t(2) * cos(nu_0);                                 % Transformation parameter
kp =  - Orbit_t(2) * sin(nu_0);                                 % Derivative of the transformation
L = [k * eye(3) zeros(3); kp * eye(3) eye(3)/(k * omega)];      % TH transformation matrix
S0 = L * S0;                                                    % TH initial boundary conditions

Phi0 = OrbitalDynamics.YA_Phi(mu, h, Orbit_t(2), 0, nu_0);      % Initial fundamental matrix
invPhi0 = (Phi0\eye(6));                                        % Inverse of the initial fundamental matrix
phi = OrbitalDynamics.YA_Phi(mu, h, Orbit_t(2), Tc, nu_f);      % Final fundamental matrix
Phi = phi * invPhi0;                                            % YA STM
r_f = Phi(1:3,:) * S0;                                          % Final dimensional position vector [m]
v_f = Phi(4:6,:) * S0;                                          % Initial dimensional velocity vector [m/s]

% Assemble the state vector
SF = [r_f; zeros(3,1)];                                         % Final conditions

%% Create the problem
% Linear problem data
params(1) = Tc;                  % TOF 
params(2) = Lc;                  % Maximum length
params(3) = Fmax;                % Maximum control authority 

params(4) = mu;                  % Gauss constant
params(5) = Orbit_t(2);          % Target orbital eccentricity
params(6) = h;                   % Angular momentum magnitude
params(7) = nu_0;                % Initial true anomaly [rad]
params(8) = nu_f;                % Final true anomaly [rad]

params(9:11) = [1 1 1].';        % Final target's docking port position in the target's body frame
params(12:14) = [1 1 1].';       % Final chaser's docking port position in the chaser's body frame // TODO: input

L = 2;                           % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 3;              % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;            % Dimension of the control vector

OptProblem = Problems.RobotLinearBerthing(S0, SF, L, StateDimension, ControlDimension, params);

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

% Dimensional space 
for i = 1:length(tau) 
    omega = mu^2 / h^3;                                            % True anomaly angular velocity
    k = 1 + Orbit_t(2) * cos(tau(i));                              % Transformation
    kp =  - Orbit_t(2) * sin(tau(i));                              % Derivative of the transformation
    L = [k * eye(3) zeros(3); kp * eye(3) eye(3)/(k * omega)];     % TH transformation matrix
    C(1:6,i) = L \ C(1:6,i);                                       % Physical space
end

% Dimensions
% C(1:3,:) = C(1:3,:) * Lc;
% tau = tau * 1000;

%% Plots
% State representation
figure
hold on
xlabel('$\nu$')
ylabel('$\mathbf{r}$')
plot(tau, C(1:3,:));
legend('$x$', '$y$', '$z$')
hold off
legend('off')
grid on;
xlim([0 tau(end)])

figure
hold on
xlabel('$\nu$')
ylabel('$\mathbf{v}$')
plot(tau, C(4:6,:));
legend('$v_x$', '$v_y$', '$v_z$')
hold off
legend('off')
grid on;
xlim([0 tau(end)])

legend('off')
grid on;
xlim([0 tau(end)])

% Propulsive acceleration plot
figure;
hold on
plot(tau, u(1:3,:), 'LineWidth', 0.3)
plot(tau, sqrt(dot(u(1:3,:),u(1:3,:),1)), 'k');
yline(Fmax, 'k--')
xlabel('$\nu$')
ylabel('$\mathbf{a}$')
legend('$a_x$', '$a_y$', '$a_z$', '$\|\mathbf{a}\|_2$', '$a_{max}$');
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