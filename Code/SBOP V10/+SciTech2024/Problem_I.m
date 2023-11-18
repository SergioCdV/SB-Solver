%% Project: SBOPT %%
% Date: 18/01/23

%% Long-range rendezvous in YA model %% 
% This script provides the solving of Problem I in URJC rendezvous laboratory,
% SciTech2024

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
baseline_flag = false;                 % Baseline (true) vs reduced solution (false)

if (baseline_flag)
    n = 20;                            % Polynomial order in the state vector expansion
    m = 100;                           % Number of sampling points
else
    n = 10;                            % Polynomial order in the state vector expansion
    m = 100;                           % Number of sampling points
end
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
% Mission constraints 
R1 = 1e3;                                          % Keep-out zone radius [m]
TOF = 16 * 3600;                                    % Maximum allowed phase time [s]

% Target orbital elements
COE = [7178e3 0.008 deg2rad(190) deg2rad(98.55) 0 0];

mu = 3.986e14;                                     % Gravitational parameter of the Earth
n = sqrt(mu/COE(1)^3);                             % Mean motion [rad/s]
K = floor(TOF/(2*pi/n));                           % Number of complete revolutions
dt = TOF - K * (2*pi/n);                           % Elapsed time in the last revolution [s]

nu_0 = OrbitalDynamics.kepler(COE);                % Initial true anomaly [rad]
COE(6) = COE(6) + n * dt;                          % Final mean anomaly [rad]
nu_f = 2*pi*K + OrbitalDynamics.kepler(COE);       % Final true anomaly [rad]

% Add linear boundary conditions
S0 = [100e3 0 75e3 -0.0542 0 -0.0418].';           % Initial conditions (relative position [m] and velocity [m/s])

% Physical parameters
Fmax = 0.5e-3;                                     % Maximum available acceleration [m/s^2]

%% Scaling 
ts = 1/n;                       % Characteristic time 
Lc = COE(1);                    % Characteristic length
Vc = Lc/ts;                     % Characteristic velocity 
gamma = Lc/ts^2;                % Characteristic acceleration 

mu = 1; 
n = 1; 

TOF = TOF / ts; 
R1 = R1 / Lc; 
Fmax = Fmax / gamma;
COE(1) = COE(1) / Lc; 
S0(1:3) = S0(1:3) / Lc; 
S0(4:6) = S0(4:6) / Vc; 

%% Final boundary conditions
% TH space transformation 
h = sqrt(mu * COE(1) * (1-COE(2)^2));                           % Target angular momentum
omega = mu^2 / h^3;                                             % True anomaly angular velocity
k = 1 + COE(2) * cos(nu_0);                                     % Transformation parameter
kp =  - COE(2) * sin(nu_0);                                     % Derivative of the transformation
L = [k * eye(3) zeros(3); kp * eye(3) eye(3)/(k * omega)];      % TH transformation matrix
S0 = L * S0;                                                    % TH initial boundary conditions

Phi0 = OrbitalDynamics.YA_Phi(mu, h, COE(2), 0, nu_0);          % Initial fundamental matrix
invPhi0 = (Phi0\eye(6));                                        % Inverse of the initial fundamental matrix
phi = OrbitalDynamics.YA_Phi(mu, h, COE(2), TOF, nu_f);         % Final fundamental matrix
Phi = phi * invPhi0;                                            % YA STM
r_f = Phi(1:3,:) * S0;                                          % Final dimensional position vector [m]
v_f = Phi(4:6,:) * S0;                                          % Initial dimensional velocity vector [m/s]

% Assemble the state vector
SF = [r_f; zeros(3,1)];                                         % Final conditions

%% Problem parameters
% Linear problem data
params(1) = TOF;                 % TOF 
params(2) = R1;                  % Maximum length
params(3) = Fmax;                % Maximum control authority 

params(4) = mu;                  % Gauss constant
params(5) = COE(2);              % Target orbital eccentricity
params(6) = h;                   % Angular momentum magnitude
params(7) = nu_0;                % Initial true anomaly [rad]
params(8) = nu_f;                % Final true anomaly [rad]

%% Create the problem
L = 2;                           % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 3;              % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;            % Dimension of the control vector

OptProblem = Problems.EllipticRendezvous(S0, SF, L, StateDimension, ControlDimension, params);

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
    k = 1 + COE(2) * cos(tau(i));                                  % Transformation
    kp =  - COE(2) * sin(tau(i));                                  % Derivative of the transformation
    L = [k * eye(3) zeros(3); kp * eye(3) eye(3)/(k * omega)];     % TH transformation matrix
    C(1:6,i) = L \ C(1:6,i);                                       % Physical space
end

%% Dimensionalization
C(1:3,:) = C(1:3,:) * Lc; 
C(4:6,:) = C(4:6,:) * Vc; 
C(7:9,:) = C(7:9,:) * gamma; 
u = u * gamma;
Fmax = Fmax * gamma;

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