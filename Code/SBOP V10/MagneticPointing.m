%% Project: SBOPT %%
% Date: 01/08/22

%% Pointing with torquers %% 
% This script provides a main interface to solve the magnetic pointing problem %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                   % Polynomial basis to be use
time_distribution = 'Legendre';       % Distribution of time intervals
n = [10 10 10 10];                    % Polynomial order in the state vector expansion
m = 100;                              % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 2;                                % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 4;                   % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;                 % Dimension of the control vector

% Problem parameters
mu = 3.986e14;                        % Gravitational parameter of the Earth

% Spacecraft orbit
a = 6900e3;                           % Orbital semimajor axis [m]
e = 0.001;                            % Orbital eccentricity
RAAN = deg2rad( 90 );                 % Mean right ascension of the ascending node  [rad]
i = deg2rad( 98 );                    % Orbital inclination     [rad]
AoP = deg2rad( 270 );                 % Mean argument of perige [rad]
theta0 = 0;                           % Reference anomaly       [rad]
COE = [a e RAAN i AoP theta0];        % Classical orbital elements

E = OrbitalDynamics.euler_matrix(COE);      % Euler matrix
D = [0 1 0; 0 0 -1; -1 0 0];                % Alignment matrix
    
% Spacecraft parameters
I = [1e-2 0 0; 0 2e-2 0; 0 0 5.5e-2];          % Inertia dyadic
T = [0.2; 0.2; 0.2];                           % Maximum moment per axis

% Operational constraints
lat = deg2rad( 40 );                           % Latitude of the point of interest
long = deg2rad( -3 );                          % Longitude of the point of interest
n_ECI = [cos(lat) * cos(long); ...
         cos(lat) * sin(long); ...
         sin(lat)];

ng = [0; 1; 0];                                % Body frame vector to which we shall point
thetanp = deg2rad( 10 );                        % Geo-pointing tolerance

n = [0; -1; 0];                                % Normal vector to the solar panels
thetasp = deg2rad( 20 );                       % Sun pointing tolerance

alpha = deg2rad( 20 );                         % Maximum axis-wise angular velocity

% Boundary conditions
q0 = [0;0;0;1];                                % Initial ECI-body frame attitude
qf = [sqrt(2)/2;0;0;sqrt(2)/2];                % Final ECI-body frame attitude
omega0 = 0.05 * rand(3,1);                     % Initial angular velocity
omegaf = zeros(3,1);                           % Final angular velocity

% Problem parameters
params = [mu COE reshape(I, 1, []) T.' alpha omega0.' omegaf.' n_ECI.' ng.' thetanp].';

% Initial and final conditions
S0 = [q0; zeros(4,1)];
SF = [qf; zeros(4,1)];

% Create the problem
OptProblem = Problems.MagneticPointing(S0, SF, L, StateDimension, ControlDimension, params);

%% Optimization
% Simple solution    
tic
[C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
toc 

% Orbital dynamics 
h = sqrt(params(1) * COE(1) * (1-COE(2)^2));            % Angular momentum 
a = COE(1) * (1-COE(2)^2) ./ (1 + COE(2)*cos(tau));     % Spacecraft position vector
dtheta = h./a.^2;                                       % Angular velocity of the anomaly

% Magnetic field
B = OrbitalDynamics.dipole_model(mu, COE, tau);

for i = 1:size(B,2)
    % Transformation to the perifocal frame 
    R = [cos(tau(i)) sin(tau(i)) 0; -sin(tau(i)) cos(tau(i)) 0; 0 0 1];
    
    % Transformation to the body frame
    S = [0 -C(3,i) C(2,i); C(3,i) 0 -C(1,i); -C(2,i) C(1,i) 0];
    A = (C(4,i)^2-C(1:3,i).'*C(1:3,i)) * eye(3) + 2 * C(1:3,i) * C(1:3,i).' - 2 * C(4,i) * S;
    B(:,i) = (A * (D * R * E).') * B(:,i);
end

% Time law
tau = cumtrapz(tau,1./dtheta);

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
xlabel('$t [s]$')
ylabel('$\mathbf{q}$')
plot(tau, C(1:4,:));
hold off
legend('off')
grid on;

figure;
ngt = repmat(ng, 1, length(tau));
for i = 1:size(C,2)
    % Transformation to the body frame
    S = [0 -C(3,i) C(2,i); C(3,i) 0 -C(1,i); -C(2,i) C(1,i) 0];
    A = (C(4,i)^2-C(1:3,i).'*C(1:3,i)) * eye(3) + 2 * C(1:3,i) * C(1:3,i).' - 2 * C(4,i) * S;

    % Pointing axis
    ngt(:,i) = A.' * ng;
end
hold on
xlabel('$t [s]$')
ylabel('$\theta_e$')
plot(tau, rad2deg( acos( abs( dot(repmat(n_ECI, 1, length(tau)), ngt, 1) ) )) );
hold off
legend('off')
grid on;

figure_orbits = figure;
omega = 2 * [C(4,:).*C(5,:)+C(3,:).*C(6,:)-C(2,:).*C(7,:)-C(1,:).*C(8,:); ...
            -C(3,:).*C(5,:)+C(4,:).*C(6,:)+C(1,:).*C(7,:)-C(2,:).*C(8,:); ...
             C(2,:).*C(5,:)-C(1,:).*C(6,:)+C(4,:).*C(7,:)-C(3,:).*C(8,:)];
hold on
xlabel('$t [s]$')
ylabel('$\mathbf{\omega}$')
plot(tau, sqrt(dot(omega,omega,1)), 'k');
yline(alpha, 'k--')
plot(tau, omega);
hold off
legend('$\|\omega\|$', '$\omega_{max}$', '$\omega_x$', '$\omega_y$', '$\omega_z$')
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
xlabel('$t [s]$')
ylabel('$\mathbf{u}$')
grid on;
