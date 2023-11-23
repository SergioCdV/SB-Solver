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
baseline_flag = true;                 % Baseline (true) vs reduced solution (false)

if (baseline_flag)
    N = 20;                            % Polynomial order in the state vector expansion
    m = 200;                           % Number of sampling points
else
    N = 10;                            % Polynomial order in the state vector expansion
    m = 200;                           % Number of sampling points
end
 
solver = Solver(basis, N, time_distribution, m);

%% Problem definition 
% Mission constraints 
R1 = 1e3;                                          % Keep-out zone radius [m]
TOF = 16 * 3600;                                   % Maximum allowed phase time [s]

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
S0 = [50e3 -10e3 7e3 -0.0542 0 -0.0418].';         % Initial conditions (relative position [m] and velocity [m/s])

% Physical parameters
Fmax = 0.5e-2;                                     % Maximum available acceleration [m/s^2]

%% Scaling 
ts = 1/n;                       % Characteristic time 
Lc = COE(1);                    % Characteristic length
Vc = Lc/ts;                     % Characteristic velocity 
gamma = Lc/ts^2;                % Characteristic acceleration 

mu = 1; 
n = 1; 

TOF = TOF / ts; 
R1 = R1 / Lc; 
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
SF = zeros(6,1);                                                % Final conditions

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
params(9) = gamma;               % Characteristic acceleration

%% Create the problem
L = 2;                           % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 3;              % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;            % Dimension of the control vector

OptProblem = Problems.EllipticRendezvous(S0, SF, L, StateDimension, ControlDimension, params);

%% Optimization
% Simple solution    
tic
[C, dV, u, t0, tf, tau, exitflag, output, P] = solver.solve(OptProblem);
toc 

% Average results 
iter = 0; 
time = zeros(1,iter);
setup.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, dV, u, t0, tf, tau, exitflag, output, P] = solver.solve(OptProblem);
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

%% Plots
% State representation
figure
subplot(1,2,1)
hold on
xlabel('$\nu$')
ylabel('$\mathbf{\rho}$ [km]', 'Interpreter','latex')
plot(tau, C(1:3,:) / 1e3);
legend('$x$', '$y$', '$z$')
hold off
grid on;
xlim([0 tau(end)])
yticklabels(strrep(yticklabels, '-', '$-$'));
xticklabels(strrep(xticklabels, '-', '$-$'));

subplot(1,2,2)
hold on
xlabel('$\nu$')
ylabel('$\dot{\mathbf{\rho}}$ [m/s]')
plot(tau, C(4:6,:) );
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$')
hold off
grid on;
xlim([0 tau(end)])
yticklabels(strrep(yticklabels, '-', '$-$'));
xticklabels(strrep(xticklabels, '-', '$-$'));

% Propulsive acceleration plot
figure;
hold on
plot(tau, u(1:3,:), 'LineWidth', 0.3)
plot(tau, sqrt(dot(u(1:3,:), u(1:3,:), 1)), 'k');
yline(Fmax, 'k--')
xlabel('$\nu$')
ylabel('$\mathbf{u}$ [m/$s^2$]')
legend('$u_r$', '$u_v$', '$u_h$', '$\|\mathbf{u}\|_2$', '$u_{max}$');
grid on;
xlim([0 tau(end)])
yticklabels(strrep(yticklabels, '-', '$-$'));
xticklabels(strrep(xticklabels, '-', '$-$'));

figure
view(3)
hold on
xlabel('$x$ [km]')
ylabel('$y$ [km]')
zlabel('$z$ [km]')
plot3(C(1,:) / 1e3, C(2,:) / 1e3, C(3,:) / 1e3);
zticklabels(strrep(zticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));
xticklabels(strrep(xticklabels, '-', '$-$'));
hold off
grid on;

%% Inertial trajectory
% Target initial orbit 
theta = linspace(tau(1), tau(end), 1e3);
r = (COE(1) * (1- COE(2)^2)) ./ (1 + COE(2) * cos(theta)) .* [cos(theta); sin(theta); zeros(1,length(theta))];
v = h/mu .* [COE(2) * sin(theta); 1 + COE(2) * cos(theta); zeros(1,length(theta))];
T = r; 
V = v;

s0 = OrbitalDynamics.coe2state(mu, COE);
h0 = cross(s0(1:3), s0(4:6));               % Angular momentum of the RSO

Q = OrbitalDynamics.euler_matrix(COE);      % Euler matrix
for i = 1:length(theta)
    T(:,i) = Q.' * r(:,i);                  % Target trajectory in inertial space

    Q1 = [cos(theta(i)) sin(theta(i)) 0; -sin(theta(i)) cos(theta(i)) 0; 0 0 1].';  % LVLH to perifocal frame
    V(:,i) = Q.' * Q1 * v(:,i);                                                     % Target hodograph in inertial space
end

% Chaser inertial orbit
S = zeros(3,length(theta));

% Rotation
o21 = -h0 / norm(h0);
for i = 1:length(theta)
    o31 = -T(:,i) / norm(T(:,i));
    o11 = cross(o21, o31);
    Q = [o11 o21 o31];

    B = PolynomialBases.Legendre().basis(N, 2 * (theta(i) / (theta(end)-theta(1))) - 1);
    S(:,i) = P * B;
    S(3,i) = S(3,i)-norm(r(:,i));
    S(:,i) = Q * S(:,i);
end

% Scaling 
S = S * Lc / 1e3;
T = T * Lc / 1e3;

figure 
view(3)
siz = 100;
hold on
scatter3(S(1,1), S(2,1), S(3,1), siz, 'b', 'Marker', 'square');
scatter3(S(1,end), S(2,end), S(3,end), siz, 'b', 'Marker', 'o');
legend('$\mathbf{s}_0$', '$\mathbf{s}_f$', 'AutoUpdate', 'off');
plot3(S(1,:), S(2,:), S(3,:), 'b'); 
plot3(T(1,:), T(2,:), T(3,:), 'r');
hold off
xlabel('$X$')
ylabel('$Y$')
zlabel('$Z$')
grid on;