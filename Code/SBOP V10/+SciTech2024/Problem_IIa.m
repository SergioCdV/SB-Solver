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
baseline_flag = true;                  % Baseline (true) vs reduced solution (false)

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
R1 = 1e3;                                          % Keep-out zone radius 1 [m]
R2 = 10;                                           % Keep-out zone radius 2 [m]
L = 1.85;                                          % Graspling reach [m]
TOF = 6 * 3600;                                    % Maximum allowed phase time [s]

Re = 6378e3;                                       % Mean Earth radius
J2 = 1.08263e-3;                                   % J2 parameter of the Earth
mu = 3.986e14;                                     % Gravitational parameter of the Earth

% Target's initial conditions (relative position [m] and velocity [m/s])
ST = [7.104981874434397e6 1.137298852087994e6 -0.756578094588272e6 -0.586250624037242e3 -1.213011751682090e3 -7.268579401702199e3].';
COE = OrbitalDynamics.state2coe(mu, ST, 'Inertial');

n = sqrt(mu/COE(1)^3);                             % Mean motion [rad/s]
K = floor(TOF/(2*pi/n));                           % Number of complete revolutions
dt = TOF - K * (2*pi/n);                           % Elapsed time in the last revolution [s]

nu_0 = mod(59.7962,2*pi);                          % Initial true anomaly [rad]
COE(6) = COE(6) + n * dt;                          % Final mean anomaly [rad]
nu_f = 2*pi*K + OrbitalDynamics.kepler(COE);       % Final true anomaly [rad]

% Initial conditions (relative position [m] and velocity [m/s])
S0 = [-365.8244397387267 -13.2607761078909 931.0422538255768 -1.4620006073131 0.6101655531506 3.2327911653314].';    

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
R2 = R2 / Lc; 
Re = Re / Lc;
L = L / Lc; 
COE(1) = COE(1) / Lc; 
S0(1:3) = S0(1:3) / Lc; 
S0(4:6) = S0(4:6) / Vc;
ST(1:3) = ST(1:3) / Lc; 
ST(4:6) = ST(4:6) / Vc; 

%% Final boundary conditions
% TH space transformation 
h = sqrt(mu * COE(1) * (1-COE(2)^2));                           % Target angular momentum

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

params(10) = R2;                 % Keep out sphere radius 2 [m]
params(11) = L;                  % Graspling reach [m]

%% Create the problem
L = 2;                           % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 3;              % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;            % Dimension of the control vector

%% Optimization (NMPC-RTI)
% Setup
options = odeset('AbsTol', 1e-22, 'RelTol', 2.25e-14);  % Integration tolerances
Ts = 60 / ts;                                           % Sampling time

% Noise 
Sigma_r = 0.1 / Lc * eye(3);                            % Relative position covariance
Sigma_v = 0.05 / Vc * eye(3);                           % Relative velocity covariance

GoOn = true;                                            % Convergence boolean
iter = 1;                                               % Initial iteration
maxIter = 1e6;                                          % Maximum number of iterations
elapsed_time = 0;                                       % Elapsed time
tof = TOF;                                              % Maximum phase time

St = zeros(maxIter,6);
C = zeros(maxIter,6);
St(1,:) = ST.';
C(1,:) = S0.';
U = [];

% Target and chaser ECI initial conditions
o21 = -cross(St(1,1:3).', St(1,4:6).');             % LVLH y+ definition
o31 = -St(1,1:3).' / norm(St(1,1:3).');             % LVLH z+ definition
o11 = cross(o21, o31);                              % LVLH x+ definition
Q = [o11 o21 o31];                                  % Rotation matrix from the LVLH to the ECI frames
y0 = [St(1,:).'; blkdiag(Q,Q) * S0 + St(1,:).'];    % ECI initial conditions

% Transformation to the TS space
omega = mu^2 / h^3;                                              % True anomaly angular velocity
k = 1 + COE(2) * cos(nu_0);                                      % Transformation parameter
kp =  - COE(2) * sin(nu_0);                                      % Derivative of the transformation
A = [k * eye(3) zeros(3); kp * eye(3) eye(3)/(k * omega)];       % TH transformation matrix
S0 = A * S0;

while (GoOn && iter < maxIter)
    % Optimization (feedback phase)
    elapsed_time = (iter-1) * Ts;
    OptProblem = Problems.EllipticRendezvous(S0, SF, L, StateDimension, ControlDimension, params);
    [~, ~, u, t0, tf, tau, exitflag, ~, P] = solver.solve(OptProblem);
    
    if (gamma * norm(u(:,1)) > Fmax)
        u(:,1) = u(:,1) / norm(u(:,1)) * Fmax/gamma;
    end
    U = [U u(:,1)];     % Control trajectory

    % Preparing phase
    solver.InitialGuessFlag = true; 
    solver.P0 = P;
    solver.maxIter = 1;
    
    % Plant dynamics 
    [~, s] = ode45(@(t,s)j2_dynamics(mu, J2, Re, t, s, U(:,end)), linspace(0, Ts, 100), y0, options);  

    % Update initial conditions and state vector
    St(iter+1,:) = s(end,1:6);                                                          % Target ECI state
    Sc = s(end,7:12);                                                                   % Chaser's ECI state
    y0 = s(end,:).';                                                                    % New initial conditions
    o21 = -cross(St(iter+1,1:3).', St(iter+1,4:6).');                                   % LVLH y+ definition
    o31 = -St(iter+1,1:3).' / norm(St(iter+1,1:3).');                                   % LVLH z+ definition
    o11 = cross(o21,o31);                                                               % LVLH x+ definition
    Q = [o11 o21 o31].';                                                                % Rotation matrix from the ECI to the LVLH frames
    rho = Sc-St(iter+1,:);                                                              % Chaser state vector in the LVLH frame
    S0 = blkdiag(Q,Q) * rho.';                                                          % Relative state vector in the LVLH frame
    C(iter+1,:) = S0.';

    % Navigation system
    S0 = mvnrnd(S0, blkdiag(Sigma_r, Sigma_v), 1).';                                    % Noisy state vector

    % Transformation to TS space
    COEo = OrbitalDynamics.state2coe(mu, St(iter+1,:).', 'Inertial');
    nu_o = OrbitalDynamics.kepler(COEo);
    h = sqrt(mu * COEo(1) * (1-COEo(2)^2));                          % Target angular momentum
    omega = mu^2 / h^3;                                              % True anomaly angular velocity
    k = 1 + COEo(2) * cos(nu_o);                                     % Transformation parameter
    kp =  - COEo(2) * sin(nu_o);                                     % Derivative of the transformation
    A = [k * eye(3) zeros(3); kp * eye(3) eye(3)/(k * omega)];       % TH transformation matrix
    S0 = A * S0;

    % Convergence 
    if (elapsed_time >= TOF || norm(C(iter+1,1:3)) <= R2)
        GoOn = false;
    else
        % Shrinking horizon
        if (elapsed_time > 0.0 * TOF)
            tof = tof - Ts; 
            params(1) = tof;

            K = floor(tof/(2*pi/n));                        % Number of complete revolutions
            dt = tof - K * (2*pi/n);                        % Elapsed time in the last revolution [s]
        end

        % New initial and final anomaly
        params(7) = nu_o;
        COEo(6) = COEo(6) + n * dt /ts;                     % Final mean anomaly [rad]
        nu_f = 2*pi*K + OrbitalDynamics.kepler(COEo);       % Final true anomaly [rad]
        params(8) = nu_f;

        % Update the number of iterations
        iter = iter + 1;
    end
end

U = [U zeros(3,1)];
t = Ts * (0:iter);      % Elapsed time vector

C = C(1:iter+1,:);
St = St(1:iter+1,:); 

%% Dimensionalization
C = C.';
C(1:3,:) = C(1:3,:) * Lc; 
C(4:6,:) = C(4:6,:) * Vc;
U = U * gamma;
t = t * ts;             

%% Plots
% State representation
figure
subplot(1,2,1)
hold on
xlabel('$t$')
ylabel('$\boldmath{\rho}$ [km]', 'Interpreter','latex')
plot(t, C(1:3,:) / 1e3);
legend('$x$', '$y$', '$z$')
hold off
grid on;
xlim([0 t(end)])

subplot(1,2,2)
hold on
xlabel('$t$')
ylabel('$\dot{\boldmath{\rho}}$ [m/s]')
plot(t, C(4:6,:) );
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$')
hold off
grid on;
xlim([0 t(end)])

% Propulsive acceleration plot
figure;
hold on
plot(t, U(1:3,:), 'LineWidth', 0.3)
plot(t, sqrt(dot(U(1:3,:), U(1:3,:), 1)), 'k');
yline(Fmax, 'k--')
xlabel('$\nu$')
ylabel('$\mathbf{u}$ [m/$s^2$]')
legend('$u_r$', '$u_v$', '$u_h$', '$\|\mathbf{u}\|_2$', '$u_{max}$');
grid on;
xlim([0 t(end)])

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

%% Auxiliary functions 
% Osculating J2 dynamics
function [dr] = j2_dynamics(mu, J2, Re, t, s, u)
    % Rotation from the LVLH to the ECI frame
    o21 = -cross(s(1:3), s(4:6));                                       % LVLH y+ definition
    o31 = -s(1:3) / norm(s(1:3));                                       % LVLH z+ definition
    o11 = cross(o21,o31);                                               % LVLH x+ definition
    Q = [o11 o21 o31];                                                  % Rotation matrix from the ECI to the LVLH frames

    % Target dynamics
    dr(1:3) = s(4:6);                                                   % Position derivative
    a = 1-1.5 * J2*(Re/norm(s(1:3)))^2 * (5*(s(3)/norm(s(1:3)))^2-1);
    b = 1-1.5 * J2*(Re/norm(s(1:3)))^2 * (5*(s(3)/norm(s(1:3)))^2-3);

    dr(4:6) = - mu * s(1:3)/norm(s(1:3))^3 .* [a; a; b];                % Velocity derivative

    % Chaser dynamics
    dr(7:9) = s(10:12);                                                 % Position derivative
    a = 1-1.5 * J2 * (Re/norm(s(7:9)))^2 * (5*(s(9)/norm(s(7:9)))^2-1);
    b = 1-1.5 * J2 * (Re/norm(s(7:9)))^2 * (5*(s(9)/norm(s(7:9)))^2-3);

    dr(10:12) = Q * u - mu * s(7:9)/norm(s(7:9))^3 .* [a; a; b];        % Velocity derivative

    dr = dr.';
end