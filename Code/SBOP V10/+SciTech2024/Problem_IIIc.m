%% Project: SBOPT %%
% Date: 18/01/23

%% Close-range rendezvous in YA model %% 
% This script provides the solving of Problem IIa) in the URJC rendezvous laboratory,
% SciTech2024

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 15;                                 % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
mu = 3.986e14;                                     % Gravitational parameter of the Earth

% Mission constraints
R2 = 10;                                           % Keep-out zone radius 2 [m]
L = 1.85;                                          % Graspling reach [m]
TOF = 6 * 3600;                                    % Maximum allowed phase time [s]

% Target's initial conditions (relative position [m] and velocity [m/s])
St = [7.104981874434397e6 1.137298852087994e6 -0.756578094588272e6 -0.586250624037242e3 -1.213011751682090e3 -7.268579401702199e3].';
COE = OrbitalDynamics.state2coe(mu, St, 'Inertial');

yt0 = [0.5 * ones(4,1); deg2rad([1 3 -1]).']; 
% ST = [-0.6952 -0.6952 -0.1294 -0.1294 0.0175 0 0];
% ST = [-0.5413 -0.4216 -0.6758 -0.2693 0.0175 -0.0079 -0.0037];
% ST = [0.2995 0.5219 -0.7845 0.1502 0.0175 -0.0158 -0.0075];
ST = [-0.3312 0.4635 -0.8219 0.0030 0.0175 0.0067 -0.0548];

n = sqrt(mu/COE(1)^3);                             % Mean motion [rad/s]
K = floor(TOF/(2*pi/n));                           % Number of complete revolutions
dt = TOF - K * (2*pi/n);                           % Elapsed time in the last revolution [s]

COEf = COE;                                        % Classical orbital elements at the final epoch
COEf(6) = COEf(6) + n * dt;                        % Final mean anomaly [rad]
nu_0 = OrbitalDynamics.kepler(COE);                % Initial true anomaly [rad]
nu_f = 2*pi*K + OrbitalDynamics.kepler(COEf);      % Final true anomaly [rad]

% Initial conditions (relative position [m] and velocity [m/s])
S0 = [467.9492284850632 -77.2962065075666 -871.9827927879848 zeros(1,3) -1.7286747525940 -0.3307280703785 5.5751101965630 zeros(1,3)].';    

% Physical parameters
Fmax = 0.5e-2;                  % Maximum available acceleration [m/s^2]
Tmax = 1e1;                     % Maximum available torque [N m]
It = diag([140 36.9 36.9]);     % Inertia tensor of the target [kg m^2]
Ic = diag([100 95.9 80]);       % Inertia tensor of the chaser [kg m^2]

theta_e = deg2rad(5);           % Error angle [rad]
theta_c = deg2rad(10);          % Cone angle [rad]
omega_max = deg2rad(100);        % Maximum chaser's angular velocity infinity norm [rad/s]
epsilon_r = deg2rad(0.1);       % Maximum relative angular velocity [rad/s]

%% Scaling 
ts = 1/n;                       % Characteristic time 
Lc = COE(1);                    % Characteristic length
Vc = Lc/ts;                     % Characteristic velocity
omega_c = n;                    % Characteristic angular velocity
gamma = Lc/ts^2;                % Characteristic acceleration 
alpha = 1/ts^2;                 % Characteristic angular acceleration

TOF = TOF / ts; 
R2 = R2 / Lc; 
L = L / Lc; 
COE([1 7])  = COE([1 7]) / Lc;
S0(1:3) = S0(1:3) / Lc; 
S0(7:9) = S0(7:9) / Vc;
Fmax = Fmax / gamma;

S0(10:12) = S0(10:12) / omega_c;
yt0(5:7) = yt0(5:7) / omega_c;
ST(5:7) = ST(5:7) / omega_c;
omega_max = omega_max / omega_c; 
epsilon_r = epsilon_r / omega_c;

Ic = Ic / Lc^2;
It = It / Lc^2;
Tmax = Tmax / (Lc * gamma);

mu = 1; 
n = 1; 

%% Final boundary conditions
% TH space transformation 
h = sqrt(mu * COE(1) * (1-COE(2)^2));      % Target angular momentum

% Assemble the state vector
SF = zeros(12,1);                          % Final conditions

%% Problem parameters
% Linear problem data
params(1) = nu_0;                % Initial true anomaly [rad]
params(2) = nu_f;                % Final true anomaly [rad]
params(3) = Fmax;                % Maximum control authority 

params(4) = mu;                  % Gauss constant
params(5) = COE(2);              % Target orbital eccentricity
params(6) = h;                   % Angular momentum magnitude

params(7) = gamma;               % Characteristic acceleration

params(9) = R2;                         % Keep out sphere radius 2 [m]
params(10) = L;                         % Graspling reach [m]

params(11:19) = reshape(Ic, 1, []);     % Inertia tensor of the chaser [kg m^2]
params(20) = Tmax;                      % Torque envelope of the chaser [N m]
params(21) = omega_max;                 % Maximum angular velocity of the chaser [rad s]
params(22:30) = reshape(It, 1, []);     % Inertia tensor of the target [kg m^2]

params(31:33) = [0 1 0];                % Boresight direction of the manipulator
params(34:36) = [0 0 1];                % Boresight direction of the graspling fixture

params(37:43) = ST;                     % Target attitude conditions

params(44) = theta_c;                   % Corridor angle restriction
params(45) = theta_e;                   % Cone angle restriction

params(46) = epsilon_r;                 % Relative angular velocity constraint 

%% Create the problem
L = 2;                                  % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 6;                     % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 6;                   % Dimension of the control vector

%% Optimization (NMPC-RTI)
% Transformation to the TS space
k = mu^2 / h^3;                                                  % True anomaly angular velocity
rho = 1 + COE(2) * cos(nu_0);                                    % Transformation parameter
drho =  - COE(2) * sin(nu_0);                                    % Derivative of the transformation
A = [rho * eye(3) zeros(3); drho * eye(3) eye(3)/(rho * k)];     % TH transformation matrix
S0([1:3 7:9]) = A * S0([1:3 7:9]);

% Definition of the problem 
OptProblem = Problems.EllipticBerthing(S0, SF, L, StateDimension, ControlDimension, params);

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
    rho = 1 + COE(2) * cos(tau(i));                                  % Transformation
    drho =  - COE(2) * sin(tau(i));                                  % Derivative of the transformation
    L = [rho * eye(3) zeros(3); drho * eye(3) eye(3)/(k * rho)];     % TH transformation matrix
    C(1:6,i) = L \ C(1:6,i);                                         % Physical space
end

% Integration of the chaser attitude 
% options = odeset('AbsTol', 1e-22, 'RelTol', 2.25e-14);
%     
% 
% [~, st] = ode45(@(t,s)euler_dynamics(It, t, s), cumtrapz(tau, Omega), yt0, options);
% st = st.';
%%
rho = 1 + COE(2) * cos(tau(1,:));    % Transformation parameter
Omega = 1 ./ (k .* rho.^2);          % Derivative of the true anomaly  
t = cumtrapz(tau, Omega);
f = (It(1,1) / It(2,2) -1);
st(5:7,:) = [yt0(5) * ones(1,length(tau)); yt0(6) * cos(f * t) + yt0(7) * sin(f * t); -yt0(6) * sin(f * t) + yt0(7) * cos(f * t)];

for i = 1:length(t)
 if (i == 1)
     st(1:4,i) = yt0(1:4);
 else
     dt = t(i)-t(i-1);
     st(1:4,i) = QuaternionAlgebra.right_isoclinic( QuaternionAlgebra.exp_map(0.5 * st(5:7,i-1) * dt) ) * st(1:4,i-1);
 end
end

%% Dimensionalization
C(1:3,:) = C(1:3,:) * Lc; 
C(7:9,:) = C(7:9,:) * Vc; 
C(13:15,:) = C(13:15,:) * gamma; 
u(1:3,:) = u(1:3,:) * gamma;
Fmax = Fmax * gamma;

st(5:7,:) = st(5:7,:) * omega_c;
epsilon_r = epsilon_r * omega_c;
omega_max = omega_max * omega_c;
Tmax = Tmax * (Lc * gamma);
u(4:6,:) = u(4:6,:) * (Lc * gamma);

Rho = C(1:3,:);             % Relative position
sigma = C(4:6,:);           % MRP
drho = C(7:9,:);            % Relative velocity
dsigma = C(10:12,:);        % Derivative of the MRP

%% Plots
% State representation
figure
subplot(1,2,1)
hold on
xlabel('$\nu$')
ylabel('$\boldmath{\rho}$ [km]', 'Interpreter','latex')
plot(tau, Rho / 1e3);
legend('$x$', '$y$', '$z$')
hold off
grid on;
xlim([tau(1) tau(end)])
yticklabels(strrep(yticklabels, '-', '$-$'));
xticklabels(strrep(xticklabels, '-', '$-$'));

subplot(1,2,2)
hold on
xlabel('$\nu$')
ylabel('$\dot{\boldmath{\rho}}$ [m/s]')
plot(tau, drho);
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$')
hold off
grid on;
xlim([tau(1) tau(end)])
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
xlim([tau(1) tau(end)])
yticklabels(strrep(yticklabels, '-', '$-$'));
xticklabels(strrep(xticklabels, '-', '$-$'));

figure
view(3)
hold on
xlabel('$x$ [km]')
ylabel('$y$ [km]')
zlabel('$z$ [km]')
plot3(Rho(1,:) / 1e3, Rho(2,:) / 1e3, Rho(3,:) / 1e3);
zticklabels(strrep(zticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));
xticklabels(strrep(xticklabels, '-', '$-$'));
hold off
grid on;

%% Attitude plots 
% Relative state representation
qr = zeros(4,length(tau));
omega_r = zeros(3,length(tau));
qc = qr;
omega = omega_r;

for i = 1:length(tau)
    q = [sigma(:,i); -1];
    B = QuaternionAlgebra.Quat2Matrix(q);
    omega(:,i) = 4 * omega_c * B.' * dsigma(:,i) / dot(q,q)^2;

    qc(:,i) = QuaternionAlgebra.MPR2Quat(1, 1, sigma(:,i), true);
    qr(:,i) = QuaternionAlgebra.right_isoclinic( qc(:,i) ) * QuaternionAlgebra.quaternion_inverse( st(1:4,i) );
    aux = QuaternionAlgebra.RotateVector( qr(:,i), st(5:7,i) );
    omega_r(:,i) = omega(:,i) - aux;
end

figure
subplot(1,2,1)
hold on
xlabel('$\nu$')
ylabel('$\mathbf{q}_r$', 'Interpreter', 'latex')
plot(tau, qr);
legend('$q_1$', '$q_2$', '$q_3$', '$q_4$')
hold off
grid on;
xlim([tau(1) tau(end)])

subplot(1,2,2)
hold on
xlabel('$\nu$')
ylabel('$\omega$ [deg/s]')
plot(tau, omega_r);
legend('$\omega_x$', '$\omega_y$', '$\omega_z$')
hold off
grid on;
xlim([tau(1) tau(end)])

figure
subplot(1,2,1)
hold on
xlabel('$\nu$')
ylabel('$\mathbf{q}_c$', 'Interpreter', 'latex')
plot(tau, qc);
legend('$q_1$', '$q_2$', '$q_3$', '$q_4$')
hold off
grid on;
xlim([tau(1) tau(end)])
yticklabels(strrep(yticklabels, '-', '$-$'));
xticklabels(strrep(xticklabels, '-', '$-$'));

subplot(1,2,2)
hold on
xlabel('$\nu$')
ylabel('$\omega$ [deg/s]')
plot(tau, deg2rad(omega));
legend('$\omega_x$', '$\omega_y$', '$\omega_z$')
hold off
grid on;
xlim([tau(1) tau(end)])
yticklabels(strrep(yticklabels, '-', '$-$'));
xticklabels(strrep(xticklabels, '-', '$-$'));

cos_theta = zeros(1,length(tau));
for i = 1:length(tau)
    b = QuaternionAlgebra.RotateVector(qr(:,i), params(34:36).');
    cos_theta(i) = dot(-params(31:33).', b);
end

figure
hold on
xlabel('$\nu$')
ylabel('$\cos\theta_e$')
plot(tau, cos_theta);
yline(cos(theta_e), 'k--');
hold off
grid on;
xlim([tau(1) tau(end)])
yticklabels(strrep(yticklabels, '-', '$-$'));
xticklabels(strrep(xticklabels, '-', '$-$'));

cos_cone = zeros(1,length(tau));
for i = 1:length(tau)
    b = QuaternionAlgebra.RotateVector(QuaternionAlgebra.quaternion_inverse( qc(:,i) ), params(31:33).');
    cos_cone(i) = dot(C(1:3,i) / norm(C(1:3,i)), b);
end

figure
hold on
xlabel('$\nu$')
ylabel('$\cos\theta_c$')
plot(tau, cos_cone);
yline(cos(theta_c), 'k--');
hold off
grid on;
xlim([tau(1) tau(end)])
yticklabels(strrep(yticklabels, '-', '$-$'));
xticklabels(strrep(xticklabels, '-', '$-$'));

% Propulsive acceleration plot
figure;
hold on
plot(tau, u(4:6,:), 'LineWidth', 0.3)
plot(tau, max(abs(u(4:6,:)), [], 1), 'k');
%yline(Tmax, 'k--')
xlabel('$\nu$')
ylabel('$\mathbf{u}$ [N m]')
legend('$\tau_x$', '$\tau_y$', '$\tau$', '$\|\mathbf{\tau}\|_\infty$');
grid on;
xlim([tau(1) tau(end)])
yticklabels(strrep(yticklabels, '-', '$-$'));
xticklabels(strrep(xticklabels, '-', '$-$'));

%% Auxiliary functions
function [dq] = euler_dynamics(I, t, s)
    c = 1e-3;

    dq(1:4) = 0.5 * QuaternionAlgebra.right_isoclinic([s(5:7); 0]) * s(1:4) + c * (1 - s(1:4).'*s(1:4)) * s(1:4);
    dq(5:7) = cross(I * s(5:7), s(5:7));

    dq = dq.';
end
