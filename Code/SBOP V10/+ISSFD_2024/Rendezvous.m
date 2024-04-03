%% Project: SBOPT %%
% Date: 31/03/24

%% Close-range rendezvous in YA model %% 
% This script provides the solving of the rendezvous problem with an
% uncooperative target for demonstration purposes in ISSFD 2024

%% Set up 
close all
clear

%% Problem definition 
% Environment definition 
mu = 3.986e14;                                     % Gravitational parameter of the Earth
Re = 6378e3;                                       % Mean Earth radius
J2 = 1.08263e-3;                                   % J2 parameter of the Earth

% Mission constraints
R1 = 1e3;                                          % Keep-out zone radius 1 [m]
R2 = 10;                                           % Keep-out zone radius 2 [m]
L = 1.85;                                          % Graspling reach [m]
Fmax = 0.5e-2;                                     % Maximum available acceleration [m/s^2]
Tmax = 0.5e-2;                                     % Maximum available torque [rad/s^2]
I = diag( [1 2 3] );                               % Inertia matrix of the chaser in the body frame [kg m2]

% Mission constraints
TOF = 6 * 3600;                                    % Maximum allowed phase time [s]

% Target's initial conditions (relative position [m], velocity [m/s], LVLH MRP, LVLH angular velocity [rad/s])
ST = [7.104981874434397e6 1.137298852087994e6 -0.756578094588272e6 -0.586250624037242e3 -1.213011751682090e3 -7.268579401702199e3].';
COE = OrbitalDynamics.ECI2COE(mu, ST, 1);

n = sqrt(mu/COE(1)^3);                             % Mean motion [rad/s]
K = floor(TOF/(2*pi/n));                           % Number of complete revolutions
dt = TOF - K * (2*pi/n);                           % Elapsed time in the last revolution [s]

COEf = COE;                                        % Classical orbital elements at the final epoch
COEf(6) = COEf(6) + n * dt;                        % Final mean anomaly [rad]

% Initial true anomaly [rad]
nu_0 = OrbitalDynamics.KeplerSolver(COE(2), COE(6)); 
if (nu_0 < 0)
    nu_0 = nu_0 + 2 * pi;
end

% Final true anomaly [rad]
nu_f = 2 * pi * K + OrbitalDynamics.KeplerSolver(COEf(2), COEf(6));  
if (nu_f < 0)
    nu_f = nu_f + 2 * pi;
end

% Initial conditions (relative position [m], velocity [m/s], LVLH MRP, LVLH angular velocity [rad/s])
S0 = [467.9492284850632 -77.2962065075666 -871.9827927879848 -1.7286747525940 -0.3307280703785 5.5751101965630].';  

%% Scaling 
ts = 1 / n;                     % Characteristic time 
Lc = COE(1);                    % Characteristic length
Vc = Lc / ts;                   % Characteristic velocity 
gamma = Lc / ts^2;              % Characteristic acceleration 
Tau = max(diag(I)) / ts^2;      % Characteristic torque

mu = 1;                         % Normalized gravitational parameter
n = 1;                          % Normalized time scale

Re = Re / Lc;                   % Normalized Earth's radius

TOF = TOF / ts;                 % Normalized mission time
R1 = R1 / Lc;                   % Normalized KOS radius
R2 = R2 / Lc;                   % Normalized KOS radius

L = L / Lc;                     % Normalized arm length

Fmax = Fmax / gamma;            % Normalized control acceleration 
Tmax = Tmax / Tau;              % Normalized torque 

% Normalized COE
COE([1 7])  = COE([1 7]) / Lc;  % Normalized target COE
COEf([1 7]) = COEf([1 7]) / Lc; % Normalized final COE

% TH space transformation 
h = sqrt(mu * COE(1) * (1-COE(2)^2));    % Target angular momentum

% Normalized initial conditions
S0(1:3) = S0(1:3) / Lc;         % Normalized initial relative position vector
S0(4:6) = S0(4:6) / Vc;         % Normalized initial relative velocity vector
ST(1:3) = ST(1:3) / Lc;         % Normalized initial target position vector
ST(4:6) = ST(4:6) / Vc;         % Normalized initial target velocity vector

%% Final boundary conditions
% Assemble the state vector
SF = zeros(6,1);               % Final reference conditions

%% Create the problem
L = 2;                           % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 3;              % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;            % Dimension of the control vector

% Linear problem data
params(1) = TOF;                 % TOF 

params(2) = mu;                  % Gauss constant
params(3) = COE(2);              % Target orbital eccentricity
params(4) = h;                   % Angular momentum magnitude
params(5) = nu_0;                % Initial true anomaly of the target [rad]
params(6) = nu_f;                % Final true anomaly of the target [rad]
params(7) = Fmax;                % Maximum control authority (linear)
% params(7) = Tmax;                % Maximum control authority (angular)
% 
% params(9) = gamma;               % Characteristic acceleration
% 
% params(10) = R1;                 % Keep out sphere radius 1 [m]
% params(11) = R2;                 % Keep out sphere radius 2 [m]
% params(12) = L;                  % Graspling reach [m]

% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
N = 10;                                % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, N, time_distribution, m);

%% Optimization (NMPC-RTI)
% Setup
options = odeset('AbsTol', 1e-22, 'RelTol', 2.25e-14);  % Integration tolerances
Ts = 200 / ts;                                          % Sampling time

% Numerial setup
GoOn = true;                                            % Convergence boolean
iter = 1;                                               % Initial iteration
maxIter = 1e3;                                          % Maximum number of iterations
elapsed_time = 0;                                       % Elapsed time

% Preallocation
St = zeros(maxIter,6);                                  % Target trajectory
C = zeros(maxIter,6);                                   % Relative trajectory
St(1,:) = ST.';                                         % Initial target state
C(1,:) = S0.';                                          % Initial relative state
U = [];                                                 % Control function

% Noise definition
Sigma_r = (0.1 / Lc)^2 * eye(3);                        % Relative position covariance
Sigma_v = (0.05 / Vc)^2 * eye(3);                       % Relative velocity covariance

% Target and chaser ECI initial conditions
Qt = OrbitalDynamics.ECI2LVLH(St(1,1:3).', St(1,4:6).', 1);      % Rotation matrix from the ECI to the LVLH frames
y0 = [St(1,:).'; St(1,:).' + blkdiag(Qt.',Qt.') * S0];           % ECI initial conditions (target - chaser)
K = 0; 

% Transformation to the TS space
omega = mu^2 / h^3;                                              % True anomaly angular velocity
k = 1 + COE(2) * cos(nu_0);                                      % Transformation parameter
kp =  - COE(2) * sin(nu_0);                                      % Derivative of the transformation
A = [k * eye(3) zeros(3); kp * eye(3) eye(3)/(k * omega)];       % TH transformation matrix
S0 = A * S0;                                                     % Initial TH relative conditions

while (GoOn && iter < maxIter)
    % Optimization (feedback phase)
    OptProblem = ISSFD_2024.RendezvousADR(S0, SF, L, StateDimension, ControlDimension, params);
    [T, ~, u, t0, tf, tau, exitflag, ~, P] = solver.solve(OptProblem);
% 
%     if (iter == 1)
%         SF = T(1:6,end);
%     end
%     
    % Control vector
%     index = gamma * sqrt(dot(u,u,1)) > Fmax;
%     u(:,index) = u(:,index) / norm(u(:,index)) * Fmax/gamma;
%     Pu = PolynomialBases.Legendre().modal_projection(u);
%     U = [U Pu * PolynomialBases.Legendre().basis(size(Pu,2)-1, 2 * (linspace(0, Ts, 100)-t0) / (tf-t0) - 1)];

    % Preparing phase
    solver.InitialGuessFlag = true; 
    solver.P0 = P;
    solver.maxIter = 1;
    
    % Plant dynamics 
    [~, s] = ode45(@(t,s)j2_dynamics(mu, J2, Re, t, s, 0, 0, 0), [0 Ts], y0, options);  
    elapsed_time = iter * Ts;

    % Update initial conditions and state vector
    y0 = s(end,:).';                                            % New initial conditions
    St(iter+1,:) = s(end,1:6);                                  % Target ECI state
    S0 = s(end,7:12)-s(end,1:6);                                % Relative state vector in the ECI frame
    C(iter+1,:) = S0;

    % Navigation system
    S0 = mvnrnd(S0.', blkdiag(Sigma_r, Sigma_v), 1).';          % Noisy state vector

    % Transformation to the TH LVLH frame 
    osc_COE = OrbitalDynamics.ECI2COE(mu, St(iter+1,:).', 1);   % Osculating target COE
    nu = OrbitalDynamics.KeplerSolver(osc_COE(2), osc_COE(6));  % Osculating true anomaly
    h = sqrt(mu * osc_COE(1) * (1 - osc_COE(2)^2));             % Osculating angular momentum

    nu_0 = OrbitalDynamics.KeplerSolver(osc_COE(2), osc_COE(6));
    if (nu_0 < 0)
        nu_0 = nu_0 + 2 * pi;
    end

    n = sqrt(mu / osc_COE(1)^3);                % Mean motion [rad/s]
    T = 2*pi / n;                               % Period of the orbit
    tf = TOF - elapsed_time;                    % Current available time
    K = max(0, floor(tf/T));                    % Number of complete revolutions
    dt = tf - K * T;                            % Elapsed time in the last revolution [s]
    osc_COEf = osc_COE;                         % Classical orbital elements at the final epoch
    osc_COEf(6) = osc_COEf(6) + max(0, n * dt); % Final mean anomaly [rad]

    nu_f = OrbitalDynamics.KeplerSolver(osc_COEf(2), osc_COEf(6));  
    if (nu_f < 0)
        nu_f = nu_f + 2 * pi;
    end

    nu_f = nu_f + 2 * pi * K;

    params(3) = osc_COE(2);      % Target orbital eccentricity
    params(4) = h;               % Angular momentum magnitude
    params(5) = nu_0;            % Initial true anomaly of the target [rad]
    params(6) = nu_f;            % Final true anomaly of the target [rad]% Update the problem parameters

    % Rotation matrix from the ECI to the LVLH frames
    Qt = OrbitalDynamics.ECI2LVLH(St(iter+1,1:3).', St(iter+1,4:6).', 1);     

    S0 = blkdiag(Qt, Qt) * S0;                                  % Initial conditions in the LVLH frame

    omega = mu^2 / h^3;                                         % True anomaly angular velocity
    k = 1 + osc_COE(2) * cos(nu);                               % Transformation parameter
    kp =  - osc_COE(2) * sin(nu);                               % Derivative of the transformation
    A = [k * eye(3) zeros(3); kp * eye(3) eye(3)/(k * omega)];  % TH transformation matrix
    S0 = A * S0;                                                % Initial conditions in the TH LVLH frame

    % Convergence 
    if (elapsed_time >= TOF)
        GoOn = false;
    else
        % Update the number of iterations
        iter = iter + 1;
    end
end

% Final processing of the results
St = St.';                   % Target trajectory
C = C.';                     % Linear relative trajectory
C = C(:,1:iter);
St = St(:,1:iter);

U = zeros(6,iter);
U = [U zeros(size(U,1),1)];  % Control vector
 
t = Ts * (0:iter-1);         % Elapsed time vector

%% Dimensionalization
t = t * ts;                     % Mission time  
C(1:3,:) = C(1:3,:) * Lc;       % Relative position 
C(4:6,:) = C(4:6,:) * Vc;       % Relative velocity
U(1:3,:) = U(1:3,:) * gamma;    % Linear control acceleration
U(4:6,:) = U(4:6,:) * Tau;      % Angular control acceleration

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
plot(linspace(0, elapsed_time, size(U,2)), U(1:3,:), 'LineWidth', 0.3)
plot(linspace(0, elapsed_time, size(U,2)), sqrt(dot(U(1:3,:), U(1:3,:), 1)), 'k');
yline(Fmax, 'k--')
xlabel('$\nu$')
ylabel('$\mathbf{u}$ [m/$s^2$]')
legend('$u_r$', '$u_v$', '$u_h$', '$\|\mathbf{u}\|_2$', '$u_{max}$');
grid on;
%xlim([0 t(end)])

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
function [dr] = j2_dynamics(mu, J2, Re, t, s, P, t0, tf)
    % LVLH transformation 
    Qt = OrbitalDynamics.ECI2LVLH(s(1:3,:), s(4:6,:), 1).';

    % Common terms
    rsquare = dot(s(1:3,:), s(1:3,:), 1);
    ReJ2 = Re^2 ./ rsquare;

    % Target dynamics
    a = 1 - 1.5 * J2 * ReJ2 .* (5 * (s(3,:).^2 ./ rsquare) - 1 );
    b = 1 - 1.5 * J2 * ReJ2 .* (5 * (s(3,:).^2 ./ rsquare) - 3 );
    dr(1:3,1) = s(4:6);                                                         % Position derivative
    dr(4:6,1) = - mu * s(1:3,:) ./ sqrt( rsquare ).^3 .* [a; a; b];             % Velocity derivative

    % Chaser dynamics
    rsquare = dot(s(7:9,:), s(7:9,:), 1);
    ReJ2 = Re^2 ./ rsquare;

    a = 1 - 1.5 * J2 * ReJ2 * ( 5 * ( s(9,:).^2 ./ rsquare) - 1 );
    b = 1 - 1.5 * J2 * ReJ2 * ( 5 * ( s(9,:).^2 ./ rsquare) - 3 );

%     tau = 2 * (t-t0) / (tf-t0) - 1;
%     B = PolynomialBases.Legendre().basis(size(P,2)-1, tau);
    u = zeros(3,1);

    dr(7:9,1) = s(10:12);                                                      % Position derivative
    dr(10:12,1) = Qt * u - mu * s(7:9,:) ./ sqrt( rsquare ).^3 .* [a; a; b];   % Velocity derivative
end