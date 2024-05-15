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
Fmax = 0.5e-2;                                     % Maximum available acceleration [m/s^2]
Tmax = 5;                                          % Maximum available torque [rad/s^2]
I = diag( [1 2 3] );                               % Inertia matrix of the chaser in the body frame [kg m2]
It = diag( [3 2 1] );                              % Inertia matrix of the target in the body frame [kg m2]
Lg = 15;                                           % Graspling distance [m]

rdoc = [1 0 0; 0 1 0].';                           % Docking ports of the chaser and target in their respective body frames [m]

vmax = 5;                                          % Maximum angular velocity of the chaser [m/s]
omega_max = 1;                                     % Maximum angular velocity of the chaser [rad/s]

% Mission constraints
TOF = 1 * 3600;                                    % Maximum allowed phase time [s]

% Target's initial conditions (relative position [m], velocity [m/s], LVLH MRP, LVLH angular velocity [rad/s])
sigma = zeros(1,3);
omega = 0 * deg2rad( [0 3.53 3.53] );
ST = [7.104981874434397e6 1.137298852087994e6 -0.756578094588272e6 -0.586250624037242e3 -1.213011751682090e3 -7.268579401702199e3 sigma omega].';

COE = OrbitalDynamics.ECI2COE(mu, ST, 1);          % Orbital elements of the target
n = sqrt(mu/COE(1)^3);                             % Mean motion [rad/s]

% Initial conditions (relative position [m], velocity [m/s], LVLH MRP, LVLH angular velocity [rad/s])
sigma = zeros(1,3);                                                             % MRP
omega = zeros(1,3);                                                             % Angular velocity in the chaser body frame
dsigma = 0.25 * omega * QuaternionAlgebra.Quat2Matrix( [sigma.'; -1] ).';       % Derivative of the MRP

S0 = 1E-1 * [467.9492284850632 -77.2962065075666 -871.9827927879848 -1.7286747525940 -0.3307280703785 5.5751101965630 sigma dsigma].';
% S0 = [1.5, 2.5 1.5, zeros(1,3) sigma, dsigma].';

%% Final boundary conditions
SF = [zeros(6,1); 1/3; 1/3; 1/3; zeros(3,1)];               % Final reference conditions

%% Scaling 
ts = 1 / n;                     % Characteristic time 
Lc = COE(1);                    % Characteristic length
Vc = Lc / ts;                   % Characteristic velocity 
gamma = Lc / ts^2;              % Characteristic acceleration 
Tau = max(diag(I)) / ts^2;      % Characteristic torque

It = It / max(diag(I));         % Normalised chaser inertia matrix
I = I / max(diag(I));           % Normalised chaser inertia matrix

TOF = TOF / ts;                 % Normalized mission time

Fmax = Fmax / gamma;            % Normalized control acceleration 
Tmax = Tmax / Tau;              % Normalized torque 

vmax = vmax / Vc;               % Normalised linear velocity 
omega_max = omega_max * ts;     % Normalised angular velocity 

Lg = Lg / Lc;                   % Normalised distance to the spacecraft

% Normalized COE
COE([1 7])  = COE([1 7]) / Lc;  % Normalized target COE

% Normalized initial conditions
S0(1:3) = S0(1:3) / Lc;         % Normalized initial relative position vector
S0(4:6) = S0(4:6) / Vc;         % Normalized initial relative velocity vector
ST(1:3) = ST(1:3) / Lc;         % Normalized initial target position vector
ST(4:6) = ST(4:6) / Vc;         % Normalized initial target velocity vector
SF(1:3) = SF(1:3) / Lc;         % Normalized reference relative position vector
SF(4:6) = SF(4:6) / Vc;         % Normalized reference relative velocity vector

S0(10:12) = S0(10:12) * ts;     % Normalized angular velocity 
ST(10:12) = ST(10:12) * ts;     % Normalized angular velocity 

SF(10:12) = SF(10:12) * ts;     % Normalized reference angular velocity 

mu = 1;                         % Normalized gravitational parameter
n = 1;                          % Normalized time scale
Re = Re / Lc;                   % Normalized Earth's radius

% TH space transformation 
h = sqrt(mu * COE(1) * (1-COE(2)^2));    % Target angular momentum

%% Create the problem
% Linear problem data
params(1) = TOF;                 % TOF 

% Initial and final anomalies
[nu_0, nu_f] = wrapp_anomaly(n, COE(2), COE(6), TOF);            

params(2) = mu;                  % Gauss constant
params(3) = COE(2);              % Target orbital eccentricity
params(4) = h;                   % Angular momentum magnitude
params(5) = nu_0;                % Initial true anomaly of the target [rad]
params(6) = nu_f;                % Final true anomaly of the target [rad]
params(7) = Fmax;                % Maximum control authority (linear)
params(8) = Tmax;                % Maximum control authority (angular)

params(9:17) = reshape(I, [], 1);       % Inertia tensor of the chaser

params(18:23) = reshape(rdoc, [], 1);   % Docking ports of the chaser and target respectively
params(24) = vmax;                      % Maximum linear velocity
params(25) = omega_max;                 % Maximum angular velocity
params(26) = Lg;                        % Graspling reach [m]

% Numerical solver definition
L = 2;                                  % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 6;                     % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 6;                   % Dimension of the control vector

basis = 'Legendre';                     % Polynomial basis to be use
time_distribution = 'Legendre';         % Distribution of time intervals
N = 7;                                  % Polynomial order in the state vector expansion
m = 30;                                 % Number of sampling points
 
solver = Solver(basis, N, time_distribution, m);

%% Optimization (NMPC-RTI)
% Setup
options = odeset('AbsTol', 1e-22, 'RelTol', 2.25e-14);  % Integration tolerances
Ts = 10 / ts;                                          % Sampling time

% Numerial setup
GoOn = true;                                            % Convergence boolean
iter = 1;                                               % Initial iteration
maxIter = ceil(TOF/Ts);                                 % Maximum number of iterations
elapsed_time = 0;                                       % Elapsed time

comp_time = zeros(1, maxIter);                          % Computational time taken for the solver

% Preallocation
St = [];                                                % Target trajectory
C = [];                                                 % Relative trajectory
U = [];                                                 % Control function
t = [];                                                 % Trajectories time

% Noise definition
Sigma_r = (0.1 / Lc)^2 * eye(3);                        % Relative position covariance
Sigma_v = (0.5 / Vc)^2 * eye(3);                        % Relative velocity covariance

% Target and chaser ECI initial conditions
Qt = OrbitalDynamics.ECI2LVLH(ST(1:3,1), ST(4:6,1), 1);                      % Rotation matrix from the ECI to the LVLH frames
y0(1:12,1) = ST(1:12,1);                                                     % Target ECI and attitude, LVLH initial conditions
y0(13:24,1) = [ST(1:6,1) + blkdiag(Qt.',Qt.') * S0(1:6); S0(7:9); omega.'];  % Chaser ECI and attitude, inertial initial conditions

omega = mu^2 / h^3;                                                          % True anomaly angular velocity
k = 1 + COE(2) * cos(nu_0);                                                  % Transformation parameter
kp =  - COE(2) * sin(nu_0);                                                  % Derivative of the transformation
y0(23) = y0(23) - (-omega * k^2);                                            % Chaser body angular veleocity with respect to the LVLH frame

% Nominal trajectory 
% [tspan, s_ref] = ode45(@(t,s)j2_dynamics(mu, J2, Re, zeros(3,2), 0, TOF, 0, 0, n, COE(2), It, I, t, s), [0 TOF], y0, options);

% Transformation to the TS space
A = [k * eye(3) zeros(3); kp * eye(3) eye(3)/(k * omega)];               % TH transformation matrix
S0(1:6) = A * S0(1:6);                                                   % Initial TH relative conditions
S0(11) = y0(23);                                                         % Body to LVLH frame angular velocity

% Preallocation of the time windows
nu_0 = [nu_0 zeros(1,maxIter-1)];
nu_f = [nu_f zeros(1,maxIter-1)];

noise = mvnrnd(zeros(1,6), blkdiag(Sigma_r, Sigma_v), maxIter).';        % Noisy state vector

while (GoOn && iter < maxIter)
    % Optimization (feedback phase)
    tic
    OptProblem = ISSFD_2024.DecoupledRendezvous(S0([1:3 7:9 4:6 10:12]), SF([1:3 7:9 4:6 10:12]), L, StateDimension, ControlDimension, params);
    [S, ~, u, t0, tf, tau, exitflag, ~, P] = solver.solve(OptProblem);

    if (exitflag ~= -2)
%         params(27:29) = S(1:3,end).';
    end

    % Control vector
    idx = dot(u(1:3,:), u(1:3,:), 1) > Fmax^2;
    u(1:3,idx) = Fmax * u(1:3,idx)./ sqrt( dot(u(1:3,idx), u(1:3,idx), 1) ); 

    for i = 1:size(tau,2)
        tau = u(4:6,i);
        if max(abs(tau)) > Tmax
            tau(abs(tau) == max(abs(tau)),:) = sign(tau(abs(tau) == max(abs(tau)),:)) * Tmax; 
            u(4:6,i) = tau;
        end
    end

    Pnew = [PolynomialBases.Legendre().modal_projection(u(1:3,:)); ...
            PolynomialBases.Legendre().modal_projection(u(4:6,:))];
    
    comp_time(iter) = toc / ts;

    % Plant dynamics 
    if ( iter == 1 )
        Pu = zeros(ControlDimension, size(u,2));
        span = linspace(0, comp_time(iter), 10);
    else
        span = linspace(t(end), t(end) + comp_time(iter), 10);
    end

    [tspan, s] = ode45(@(t,s)j2_dynamics(mu, J2, Re, Pu, t0, tf, Fmax, Tmax, n, params(3), It, I, t-elapsed_time, s), span, y0, options); 

    % Shadow transformation
    idx = dot(s(:,7:9), s(:,7:9), 2) > 1;
    s(idx,7:9) = -s(idx,7:9) ./ dot(s(idx,7:9), s(idx,7:9), 2);

    idx = dot(s(:,19:21), s(:,19:21), 2) > 1;
    s(idx,19:21) = -s(idx,19:21) ./ dot(s(idx,19:21), s(idx,19:21), 2);

    if (iter == 1)
        t = tspan;

        % Save the trajectory
        St = [St; s(:,1:12)];                                     % Target ECI and attitude state
        C = [C; [s(:,13:18)-s(:,1:6) s(:,19:24)] ];               % Output trajectory

    else
        t = [t; tspan(2:end)];

        % Save the trajectory
        St = [St; s(2:end,1:12)];                                 % Target ECI and attitude state
        C = [C; [s(2:end,13:18)-s(2:end,1:6) s(2:end,19:24)] ];   % Output trajectory
    end

    % Control law
    tau = zeros(1, length(tspan));
    for i = 1:size(tspan,1)
        [~, nu] = wrapp_anomaly(n, params(3), t0, tspan(i)-elapsed_time);
        if (tf <= t0)
            tau(i) = -1;
        else
            tau(i) = 2 * (nu-t0) / (tf-t0) - 1;                    % Evaluation point for the controller
        end
        tau(i) = max(-1, min(1, tau(i)));
    end

    u_aux = Pu * PolynomialBases.Legendre().basis( size(Pu,2)-1, tau );

    idx = dot(u_aux(1:3,:), u_aux(1:3,:), 1) > Fmax^2;
    u_aux(1:3,idx) = Fmax * u_aux(1:3,idx) ./ sqrt( dot(u_aux(1:3,idx), u_aux(1:3,idx), 1) ); 

    for i = 1:size(tau,2)
        tau = u_aux(4:6,i);
        if max(abs(tau)) > Tmax
            tau(abs(tau) == max(abs(tau)),:) = sign(tau(abs(tau) == max(abs(tau)),:)) * Tmax; 
            u_aux(4:6,i) = tau;
        end
    end

    if (iter == 1)
        U = u_aux;
    else
        U = [ U u_aux(:,2:end) ];
    end

    % New controller
    if (exitflag ~= -2)
        Pu = Pnew;
    end
    
    % Plant dynamics 
    y0 = s(end,:).';
    [tspan, s] = ode45(@(t,s)j2_dynamics(mu, J2, Re, Pu, t0, tf, Fmax, Tmax, n, params(3), It, I, t-elapsed_time, s), linspace(t(end), t(end) + Ts, 10), y0, options);  

    % Control law
    tau = zeros(1, length(tspan));
    for i = 1:size(tspan,1)
        [~, nu] = wrapp_anomaly(n, params(3), t0, tspan(i)-elapsed_time);
        if (tf <= t0)
            tau(i) = -1;
        else
            tau(i) = 2 * (nu-t0) / (tf-t0) - 1;                                % Evaluation point for the controller
        end
        tau(i) = max(-1, min(1, tau(i)));
    end

    u_aux = Pu * PolynomialBases.Legendre().basis( size(Pu,2)-1, tau );

    idx = dot(u_aux(1:3,:), u_aux(1:3,:), 1) > Fmax^2;
    u_aux(1:3,idx) = Fmax * u_aux(1:3,idx) ./ sqrt( dot(u_aux(1:3,idx), u_aux(1:3,idx), 1) ); 

    for i = 1:size(tau,2)
        tau = u_aux(4:6,i);
        if max(abs(tau)) > Tmax
            tau(abs(tau) == max(abs(tau)),:) = sign(tau(abs(tau) == max(abs(tau)),:)) * Tmax; 
            u_aux(4:6,i) = tau;
        end
    end

    U = [ U u_aux(:,2:end) ];

    % Shadow transformation
    idx = dot(s(:,7:9), s(:,7:9), 2) > 1;
    s(idx,7:9) = -s(idx,7:9) ./ dot(s(idx,7:9), s(idx,7:9), 2);

    idx = dot(s(:,19:21), s(:,19:21), 2) > 1;
    s(idx,19:21) = -s(idx,19:21) ./ dot(s(idx,19:21), s(idx,19:21), 2);

    % Save the trajectory
    t = [t; tspan(2:end)];
    St = [St; s(2:end,1:12)];                                 % Target ECI and attitude state
    C = [C; [s(2:end,13:18)-s(2:end,1:6) s(2:end,19:24)] ];   % Output trajectory

    elapsed_time = t(end);

    % Update initial conditions and state vector
    y0 = s(end,:).';                                      % New initial conditions
    S0(1:6) = y0(13:18,1) - y0(1:6,1);                    % Relative state vector in the ECI frame
    S0(7:12) = y0(19:24,1);                               % Chaser attitude state

    % Navigation system
%     S0(1:6,1) = S0(1:6,1) + noise(:,iter);

    % Transformation to the TH LVLH frame 
    osc_COE = OrbitalDynamics.ECI2COE(mu, y0, 1);               % Osculating target COE
    n = sqrt(mu / osc_COE(1)^3);                                % Osculating mean motion
    h = sqrt(mu * osc_COE(1) * (1 - osc_COE(2)^2));             % Osculating angular momentum

    [nu_0(iter+1), nu_f(iter+1)] = wrapp_anomaly(n, osc_COE(2), osc_COE(6), max(0, TOF-elapsed_time));

    nu = nu_0(iter+1);                                          % Osculating true anomaly
    params(3) = osc_COE(2);                                     % Target orbital eccentricity
    params(4) = h;                                              % Angular momentum magnitude
    params(5) = nu_0(iter+1);                                   % Initial true anomaly of the target [rad]
    params(6) = nu_f(iter+1);                                   % Final true anomaly of the target [rad]% Update the problem parameters

    % Rotation matrix from the ECI to the LVLH frames
    Qt = OrbitalDynamics.ECI2LVLH(y0(1:3,1), y0(4:6,1), 1);     

    S0(1:6) = blkdiag(Qt, Qt) * S0(1:6);                        % Initial conditions in the LVLH frame

    omega = mu^2 / h^3;                                         % True anomaly angular velocity
    k = 1 + osc_COE(2) * cos(nu);                               % Transformation parameter
    kp =  - osc_COE(2) * sin(nu);                               % Derivative of the transformation
    A = [k * eye(3) zeros(3); kp * eye(3) eye(3)/(k * omega)];  % TH transformation matrix
    S0(1:6) = A * S0(1:6);                                      % Initial conditions in the TH LVLH frame

    S0(10:12) = 0.25 * QuaternionAlgebra.Quat2Matrix( [S0(7:9,1); -1] ) * S0(10:12);       % Derivative of the MRP

    % Convergence 
    if (elapsed_time >= TOF)
        GoOn = false;
    else
        % Update the number of iterations
        iter = iter + 1;
    end
end

% Final processing of the results
St = St.';                      % Target trajectory
C = C.';                        % Linear relative trajectory and chaser attitude

nu_0 = nu_0(1:iter);
nu_f = nu_f(1:iter);
 
%% Dimensionalization
t = t * ts;                     % Mission time  
C(1:3,:) = C(1:3,:) * Lc;       % Relative position 
C(4:6,:) = C(4:6,:) * Vc;       % Relative velocity
C(10:12,:) = C(10:12,:) / ts;   % Normalised derivative of the MRP 

U(1:3,:) = U(1:3,:) * gamma;    % Linear control acceleration
U(4:6,:) = U(4:6,:) * Tau;      % Angular control acceleration

comp_time = comp_time * ts;

%% Save results 
save Rendezvous_Uncoupled_ISSFD2024.mat;

Trajectory = [t / ts C(1:3,:).' / Lc zeros(size(C,2),3) ones(size(C,2),1)];
csvwrite('RVD_example_I.csv', Trajectory);

%% Plots
% State representation
figure
subplot(1,2,1)
hold on
xlabel('$t$ [s]')
ylabel('$\boldmath{\rho}$ [m]', 'Interpreter', 'latex')
plot(t, C(1:3,:));
yline(Lg * Lc, 'k--')
legend('$x$', '$y$', '$z$', '$L_g$', 'Autoupdate', 'off')
yline(-Lg * Lc, 'k--')
hold off
grid on;
xlim([0 t(end)])
yticklabels(strrep(yticklabels, '-', '$-$'));

subplot(1,2,2)
hold on
xlabel('$t$')
ylabel('$\dot{\boldmath{\rho}}$ [m/s]')
plot(t, C(4:6,:) );
yline(vmax * Vc, 'k--')
legend('$\dot{x}$', '$\dot{y}$', '$\dot{z}$', '$v_{max}$')
hold off
grid on;
xlim([0 t(end)])

figure
subplot(1,2,1)
hold on
xlabel('$t$ [s]')
ylabel('$\boldmath{\sigma}$', 'Interpreter','latex')
plot(t, C(7:9,:));
yline(SF(7:9), 'k--')
legend('$\sigma_x$', '$\sigma_y$', '$\sigma_z$')
hold off
grid on;
xlim([0 t(end)])
yticklabels(strrep(yticklabels, '-', '$-$'));

subplot(1,2,2)
hold on
xlabel('$t$')
ylabel('$\omega$ [rad/s]')
plot(t, C(10:12,:) );
legend('$\omega_x$', '$\omega_y$', '$\omega_z$')
hold off
grid on;
xlim([0 t(end)])

% Propulsive acceleration plot
figure;
hold on
plot(t, U(1:3,:), 'LineWidth', 0.3)
plot(t, sqrt(dot(U(1:3,:), U(1:3,:), 1)), 'k');
yline(Fmax * gamma, 'k--')
xlabel('$t$ [s]')
ylabel('$\mathbf{u}$ [m/$s^2$]')
legend('$u_r$', '$u_v$', '$u_h$', '$\|\mathbf{u}\|_2$', '$u_{max}$');
grid on;
xlim([0 t(end)])

% Propulsive acceleration plot
figure;
hold on
plot(t, U(4:6,:), 'LineWidth', 0.3)
plot(t, max(abs(U(4:6,:))), 'k');
yline(Tmax * Tau, 'k--')
xlabel('$\nu$')
ylabel('$\mathbf{\tau}$ [Nm]')
legend('$\tau_x$', '$\tau_y$', '$\tau_z$', '$\|\mathbf{\tau}\|_{\infty}$', '$\tau_{max}$');
grid on;
xlim([0 t(end)])

figure
view(3)
hold on
xlabel('$x$ [km]')
ylabel('$y$ [km]')
zlabel('$z$ [km]')
plot3(C(1,:) / 1e3, C(2,:) / 1e3, C(3,:) / 1e3);
% zticklabels(strrep(zticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));
% xticklabels(strrep(xticklabels, '-', '$-$'));
hold off
grid on;

if false
    r1 = [cos(nu_0); sin(nu_0)];
    r2 = [cos(nu_f); sin(nu_f)];
    figure 
    grid on;
    xlim([-1 1])
    ylim([-1 1])
    for i = 1:size(nu_0,2)
        hold on
        J = plot(r1(1,i), r1(2,i), 'or');
        H = plot(r2(1,i), r2(2,i), 'ob');
        drawnow;
        pause(2)
        delete(J);
        delete(H)
    end
end

figure 
stem(1:iter, comp_time(1:iter) );
grid on
xlabel('Comp. time [s]')
ylabel('Time step')

%% Auxiliary functions 
% CW dynamics 
function [drho] = CW_dynamics(t,s,P,t0,tf,Fmax)

    drho(1:3,:) = s(4:6,:);
    tau = 2 * (t) / (tf-t0) - 1;                                % Evaluation point for the controller
    tau = min(1, max(-1, tau));

    B = PolynomialBases.Legendre().basis(size(P,2)-1, tau);
    u = P * B;

    if norm(u) > Fmax
        u = Fmax * u / norm(u); 
    end

    drho(4:6,:) = u + [2 * s(6,:); -s(2,:); 3 * s(3,:) - 2 * s(4,:)];

end

% Osculating J2 dynamics
function [dr] = j2_dynamics(mu, J2, Re, P, t0, tf, Fmax, Tmax, n, e, It, Ic, t, s)
    % LVLH transformation 
    Qt = OrbitalDynamics.ECI2LVLH(s(1:3,:), s(4:6,:), 1);

    % Common terms
    rsquare = dot(s(1:3,:), s(1:3,:), 1);
    ReJ2 = Re^2 ./ rsquare;

    % Target linear dynamics
    a = 1 - 1.5 * J2 * ReJ2 .* (5 * (s(3,:).^2 ./ rsquare) - 1 );
    b = 1 - 1.5 * J2 * ReJ2 .* (5 * (s(3,:).^2 ./ rsquare) - 3 );
    dr(1:3,1) = s(4:6);                                                         % Position derivative
    dr(4:6,1) = - mu * s(1:3,:) ./ sqrt( rsquare ).^3 .* [a; a; b];             % Velocity derivative

    % Target attitude dynamics 
    sigma_t = s(7:9,1);
    if ( dot(sigma_t, sigma_t, 1) > 1 )
        sigma_t = - sigma_t / dot(sigma_t, sigma_t, 1);
    end

    sigma_t = [sigma_t; -1];

    dr(7:9,1) = 0.25 * QuaternionAlgebra.Quat2Matrix( sigma_t ) * s(10:12,1);   % Target attitude kinematics

    h = cross(s(1:3,:), s(4:6,:));
    dtheta = sqrt(dot(h,h,1)) ./ rsquare;
    omega_t =  s(10:12,1) + [0; -dtheta; 0];
    vt = Qt * s(4:6,:);
    alpha = -2 * dtheta .* (-vt(3,:) ./ vt(1,:)) * [0; -1; 0];

    qm = QuaternionAlgebra.MPR2Quat(1, 1, s(7:9,:), true);
    rb = QuaternionAlgebra.RotateVector(qm, Qt * s(1:3,:));

    dr(10:12,1) = It \ (-It * alpha + cross(It * omega_t, omega_t) + mu ./ sqrt( rsquare ).^5 .* cross(rb, Ic * rb));  % Target attitude dynamics

    % Chaser linear dynamics
    rsquare = dot(s(13:15,:), s(13:15,:), 1);
    ReJ2 = Re^2 ./ rsquare;

    a = 1 - 1.5 * J2 * ReJ2 * ( 5 * ( s(15,:).^2 ./ rsquare) - 1 );
    b = 1 - 1.5 * J2 * ReJ2 * ( 5 * ( s(15,:).^2 ./ rsquare) - 3 );

    % Controller
    [~, nu] = wrapp_anomaly(n, e, t0, t);

    if (tf <= t0)
        tau = -1;
    else
        tau = 2 * (nu-t0) / (tf-t0) - 1;                                % Evaluation point for the controller
    end
    tau = max(-1, min(1, tau));

    B = PolynomialBases.Legendre().basis(size(P,2)-1, tau);
    u = P * B;

    % Linear control
    if norm(u(1:3,:)) > Fmax
        u(1:3,:) = Fmax * u(1:3,:) / norm(u(1:3,:)); 
    end

    qm = QuaternionAlgebra.MPR2Quat(1, 1, s(19:21,:), true);
    inv_qm = QuaternionAlgebra.quaternion_inverse( qm );
    thrust = u(1:3,:);%QuaternionAlgebra.RotateVector(inv_qm, 0 * u(1:3,:));

    dr(13:15,1) = s(16:18);                                                             % Position derivative
    dr(16:18,1) = Qt.' * thrust - mu * s(13:15,:) ./ sqrt( rsquare ).^3 .* [a; a; b];   % Velocity derivative

    % Torque
    tau = u(4:6,:);
    if max(abs(tau)) > Tmax
        tau(abs(tau) == max(abs(tau)),:) = sign(tau(abs(tau) == max(abs(tau)),:)) * Tmax; 
    end

    % Chaser attitude dynamics
    sigma_c = s(19:21,1);
    if ( dot(sigma_c, sigma_c, 1) > 1 )
        sigma_c = - sigma_c / dot(sigma_c, sigma_c, 1);
    end

    sigma_c = [sigma_c; -1];

    omega_t = s(22:24,1) + [0; -dtheta; 0];
    dr(19:21,1) = 0.25 * QuaternionAlgebra.Quat2Matrix( sigma_c ) * s(22:24,1);           % Target attitude kinematics

    rb = QuaternionAlgebra.RotateVector(qm, Qt * s(13:15,:));

    dr(22:24,1) = Ic \ (tau - Ic * alpha + cross(Ic * omega_t, omega_t) + mu ./ sqrt( rsquare )^5 .* cross(rb, Ic * rb));   % Target attitude dynamics
end

% Attitude dynamics 
function [dq] = attitude_dynamics(I, t, s)
    % Variables 
    q = s(1:4,1);       % Attitude quaternion 
    omega = s(5:7,1);   % Angular velocity 

    % Kinematics 
    dq(1:4,1) = 0.5 * QuaternionAlgebra.right_isoclinic([omega; 0]) * q + 1E-3 * (1 - dot(q,q,1)) * q;

    % Dynamics 
    dq(5:7,1) = I \ cross(I * omega, omega);
end

% Compute the final true anomaly considering multiple revolutions 
function [nu_0, nu_f] = wrapp_anomaly(n, e, M, t)
    % Initial true anomaly [rad]
    nu_0 = OrbitalDynamics.KeplerSolver(e, M); 

    % Final true anomaly
    T = 2*pi/ n;                         % Orbital period
    K = floor(t / T);                    % Number of complete revolutions
    dt = t - K * (2*pi/n);               % Elapsed time in the last revolution [s]
    Mf = M + n * dt;                     % Final mean anomaly [rad]
    
    % Final true anomaly [rad]
    nu_f = OrbitalDynamics.KeplerSolver(e, Mf);  
    
    dnu = nu_f - nu_0; 
    if (dnu < 0)
        dnu = 2 * pi + dnu;
    end
    
    nu_f = 2 * pi * K + nu_0 + dnu;
end