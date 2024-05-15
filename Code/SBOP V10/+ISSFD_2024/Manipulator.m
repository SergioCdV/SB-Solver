%% Project: SBOPT %%
% Date: 31/03/24

%% Capture phase %% 
% This script provides the solving of the manipulator problem for
% demonstration purposes in ISSFD 2024

%% Set up 
close all
clear

%% Problem definition 
% Mission constraints
Tc = 300;                              % Characteristic time [s]
T0 = 0;                                % Initial clcok
TOF = 0.3 * Tc;                          % Mission TOF
Omega_max = [pi; 2*pi];                % Maximum angular velocity of the joints [rad/s]
vmax = 2;                              % Maximum linear velocity [m]
omega_max = deg2rad(20);               % Maximum angular velocity of the end-effector [rad/s]
epsilon = 1e-4^2;                      % Numerical tolerance for the Jacobian determinant

% Linear problem data
params(1) = T0;                        % Initial time [s]
params(2) = TOF;                       % Final time [s]
params(3:4) = Omega_max;               % Maximum control authority 
params(5) = vmax;                      % Constrain on the linear velocity
params(6) = omega_max;                 % Constraint on the angular velocity
params(7) = epsilon;                   % Numerical tolerance for the Jacobian determinant

params(8:10) =  [0; 0; 0];             % Base of the DH parameters
params(11:16) = [0 0 0 0 0 0].';
params(17:22) = [pi/2 0 0 pi/2 -pi/2 0].';
params(23:28) = [0 0 0 0 0 0].';
params(29:34) = [0 -0.24365 -0.21325 0 0 0].';
params(35:40) = [0.15185 0 0 0.1124 0.08535 0.0921].';
params(41:46) = ones(6,1);

% Final conditions of the end-effector (reference position and attitude)
sigma = [0.33;0.01;0.33];
params(47:58) = [0.2; -0.05; 0.01; sigma; zeros(3,1); deg2rad([0 0 5]).'].';

S0 = [0 deg2rad(-151.31) deg2rad(142.22) deg2rad(-145) deg2rad(89.37) deg2rad(125)].';
%S0 = [0 deg2rad(-135) pi/2 deg2rad(-135) pi/2 pi/2].';
SF = [0 -3*pi/4 +pi/2 -3*pi/4 pi/2 0].';

%% Numerical solver definition 
L = 1;                                 % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 6;                    % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 6;                  % Dimension of the control vector

basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 5;                                % Polynomial order in the state vector expansion
m = 20;                                % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Optimization (NMPC-RTI)
% Setup
options = odeset('AbsTol', 1e-22, 'RelTol', 2.25e-14);  % Integration tolerances
Ts = 2;                                                 % Sampling time [s]

% Numerical setup 
GoOn = true;                                            % Convergence boolean
iter = 1;                                               % Initial iteration
maxIter = ceil(TOF/Ts) + 1;                             % Maximum number of iterations
elapsed_time = 0;                                       % Elapsed time

comp_time = zeros(1, maxIter);

% Preallocation
C = [];                                                 % Relative trajectory
U = [];                                                 % Control function
t = [];                                                 % Trajectories time

% Noise 
Sigma = deg2rad(0.1)^2 * eye(6);                        % Joint state covariance [rad]

% Initial conditions
y0 = S0.';

% DH parameters of the robot
DH_parameters.base =   reshape(params(8:10), 1, 3).'; 
DH_parameters.theta =  reshape(params(11:16), 1, StateDimension).';
DH_parameters.alpha =  reshape(params(17:22), 1, StateDimension).';
DH_parameters.offset = reshape(params(23:28), 1, StateDimension).';
DH_parameters.a =      reshape(params(29:34), 1, StateDimension).';
DH_parameters.d =      reshape(params(35:40), 1, StateDimension).';
DH_parameters.type =   reshape(params(41:46), 1, StateDimension).';

while (GoOn && iter < maxIter)
    % Optimization (feedback phase)
    tic
    OptProblem = ISSFD_2024.RobotKinematics(S0, SF, L, StateDimension, ControlDimension, params);
    [S, dV, u, t0, tf, tau, exitflag, ~, P] = solver.solve(OptProblem);

    % Control vector
%     idx2 = abs(u(1:4,:)) == max(abs(u(1:4,:)), [], 1);
%     idx = u(idx2) > Omega_max(1);
%     u(idx2(idx)) = repmat(Omega_max(1), sum(idx2(idx),1), sum(idx));
%     
%     idx2 = abs(u(5:6,:)) == max(abs(u(5:6,:)), [], 1);
%     idx = u(idx2) > Omega_max(2);
%     u(idx2(idx)) = repmat(Omega_max(2), sum(idx2,1), sum(idx));
    
    Pnew = PolynomialBases.Legendre().modal_projection(u);

    comp_time(iter) = toc;
    
    % Plant dynamics
    if ( iter == 1 )
        Pu = zeros(ControlDimension, size(u,2));
        span = linspace(0, comp_time(iter), 10);
    else
        span = linspace(t(end), t(end) + comp_time(iter), 10);
    end

    [tspan, s] = ode45(@(t,s)robot_kinematics(t, s, Pu, t0, tf, Omega_max), span, y0, options);
    
    % Control law
    tau = zeros(1, length(tspan));
    for i = 1:size(tspan,1)
        if (tf <= t0)
            tau(i) = -1;
        else
            tau(i) = 2 * (tspan(i)-t0) / (tf-t0) - 1;     % Evaluation point for the controller
        end
        tau(i) = min(1, max(-1, tau(i)));
    end

    u_aux = Pu * PolynomialBases.Legendre().basis( size(Pu,2)-1, tau );

    ds = u_aux.';

    if (iter == 1)
        t = tspan;

        % Save the trajectory
        C = [C; [s ds] ];            % Output trajectory

    else
        t = [t; tspan(2:end)];

        % Save the trajectory
        C = [C; [s(2:end,:) ds(2:end,:)] ];   % Output trajectory
    end

%     idx2 = abs(u_aux(1:4,:)) == max(abs(u_aux(1:4,:)), [], 1);
%     idx = u_aux(idx2) > Omega_max(1);
%     u_aux(idx2(idx)) = repmat(Omega_max(1), sum(idx2,1), 1);
%     
%     idx2 = abs(u_aux(5:6,:)) == max(abs(u_aux(5:6,:)), [], 1);
%     idx = u_aux(idx2) > Omega_max(2);
%     u_aux(idx2(idx)) = repmat(Omega_max(2), sum(idx2,1), sum(idx));

    if (iter == 1)
        U = u_aux;
    else
        U = [ U u_aux(:,2:end) ];
    end

    % New controller
    if (exitflag ~= -2)
        Pu = Pnew;

        if (iter == 1)
            params(59:64) = mod(S(1:6,end).', 2*pi);
        end
    
    end

    % Plant dynamics
    y0 = mod(s(end,:), 2*pi);
    [tspan, s] = ode45(@(t,s)robot_kinematics(t, s, Pu, t0, tf, Omega_max), linspace(t(end), t(end) + Ts, 10), y0, options); 

    tau = zeros(1, length(tspan));
    for i = 1:size(tspan,1)
        if (tf <= t0)
            tau(i) = -1;
        else
            tau(i) = 2 * (tspan(i)-t0) / (tf-t0) - 1;     % Evaluation point for the controller
        end
        tau(i) = min(1, max(-1, tau(i)));
    end

    u_aux = Pu * PolynomialBases.Legendre().basis( size(Pu,2)-1, tau );
%     idx2 = abs(u_aux(1:4,:)) == max(abs(u_aux(1:4,:)), [], 1);
%     idx = u_aux(idx2) > Omega_max(1);
%     u_aux(idx2(idx)) = repmat(Omega_max(1), sum(idx2,1), 1);
%     
%     idx2 = abs(u_aux(5:6,:)) == max(abs(u_aux(5:6,:)), [], 1);
%     idx = u_aux(idx2) > Omega_max(2);
%     u_aux(idx2(idx)) = repmat(Omega_max(2), sum(idx2,1), sum(idx));

    ds = u_aux.';

    U = [ U u_aux(:,2:end) ];

    % Save the trajectory
    t = [t; tspan(2:end)];
    C = [C; [s(2:end,:) ds(2:end,:)] ];               % Output trajectory

    elapsed_time = t(end);

    % Update initial conditions
    y0 = mod(s(end,:), 2*pi);
    S0 = y0.';
%     params(1) = elapsed_time;                        % Initial time [s]

    % Navigation system
    noise = mvnrnd(zeros(1,size(y0,2)), Sigma, 1).';          % Noisy state vector
    S0 = S0 + noise;

    % Convergence 
    if (elapsed_time >= TOF)
        GoOn = false;
    else
        % Update the number of iterations
        iter = iter + 1;
    end
end

% Final processing of the results
C = C.';                        % Output trajectory

%% Analysis
% Check for singularities 
detJ = zeros(1, length(t)); 
v = zeros(StateDimension, length(t));
r_eff = zeros(3, length(t));
sigma = zeros(3, length(t));

for i = 1:length(t)
    % Compute the Jacobian 
    [T, J] = ISSFD_2024.RobotKinematics.Kinematics(OptProblem.StateDim, @(j,s)ISSFD_2024.RobotKinematics.ur3_dkinematics(DH_parameters, j, s), C(:,i));

    % Compute the determinant of the Jacobian
    detJ(i) = det(J);

    % Compute the frame velocities
    v(:,i) = J * C(7:12,i);

    % End-effector system
    r_eff(:,i) = T(1:3,end);                     % End effector position 

    R = T(1:3,end-2:end);                   % Final rotation matrix
    q = QuaternionAlgebra.Matrix2Quat(R);   % Associated quaternion
    
    sigma(:,i) = QuaternionAlgebra.MPR2Quat(1, 1, q, false);
end

% Final position and attitude error
r_ref = params(47:49).';                % Reference position
sigma_ref = params(50:52).';            % Reference attitude
v_ref = params(53:58).';                % Reference velocity

res(1:3,1) = r_eff(:,end) - r_ref;      % Error to the reference position
res(4:6,1) = sigma(:,end) - sigma_ref;  % Error to the reference attitude 
res(7:12,1) = v(:,end) - v_ref;         % Error to the velocities

% Close loop cost 
close_cost = trapz(linspace(0, elapsed_time, size(U,2)), dot(U,U,1));

C(1:6,:) = mod(C(1:6,:), 2*pi);

%% Save results 
save Manipulator_ISSFD2024.mat;

%% Plots
% State representation
figure
hold on
xlabel('$t$ [s]')
ylabel('$\mathbf{\theta}$ [deg]')
plot(t, rad2deg(C(1:6,:)));
legend('$\theta_1$', '$\theta_2$', '$\theta_3$', '$\theta_4$', '$\theta_5$', '$\theta_6$')
hold off
grid on;
xlim([0 t(end)])
yticklabels(strrep(yticklabels, '-', '$-$'));

figure
hold on
xlabel('$t$ [s]')
ylabel('$\mathbf{r}_e$ [m]')
plot(t, r_eff);
yline(params(47:49), 'k--')
legend('$r_x$', '$r_y$', '$r_z$', '$r_{ref}^x$', '$r_{ref}^y$', '$r_{ref}^z$')
hold off
grid on;
xlim([0 t(end)])
yticklabels(strrep(yticklabels, '-', '$-$'));

figure
hold on
xlabel('$t$ [s]')
ylabel('$\mathbf{r}_e$ [m]')
plot(t, sigma);
yline(params(50:52), 'k--')
legend('$\sigma_x$', '$\sigma_y$', '$\sigma_z$', '$\sigma_{ref}^x$', '$\sigma_{ref}^y$', '$\sigma_{ref}^z$')
hold off
grid on;
xlim([0 t(end)])
yticklabels(strrep(yticklabels, '-', '$-$'));

figure
hold on
xlabel('$t$ [s]')
ylabel('$\mathbf{r}_e$ [m]')
plot(t, v(4:6,:));
yline(params(56:58), 'k--')
legend('$\omega_x$', '$\omega_y$', '$\omega_z$', '$\omega_{ref}^x$', '$\omega_{ref}^y$', '$\omega_{ref}^z$')
hold off
grid on;
xlim([0 t(end)])
yticklabels(strrep(yticklabels, '-', '$-$'));

% Propulsive acceleration plot
figure;
hold on
plot(t, max(abs(U(1:3,:)), [], 1), 'r');
plot(t, max(abs(U(4:6,:)), [], 1), 'b');
yline(Omega_max(1), 'r--')
yline(Omega_max(2), 'b--')
xlabel('$t$')
ylabel('$\mathbf{u}$ [rad/s]')
legend('$\|\omega^{1:3}\|_{\infty}$', '$\|\omega^{4:6}\|_{\infty}$', '$\omega_{max,1}$', '$\omega_{max,2}$');
grid on;
xlim([0 t(end)])

% Frame velocities plot
figure;
hold on
subplot(2,1,1)
plot(t, max(abs(v(4:6,:)), [], 1), 'b');
yline(params(6), 'r--')
legend('$\|\mathbf{\omega}\|_{\infty}$', '$\omega_{max}$');
grid on;
subplot(2,1,2)
plot(t, max(abs(v(1:3,:)), [], 1), 'r');
yline(params(5), 'b--')
xlabel('$t$')
legend('$\|\dot{\mathbf{d}}\|_{\infty}$', '$v_{max}$');
grid on;
xlim([0 t(end)])

figure
hold on
xlabel('$t$')
ylabel('det $J$')
plot(t, detJ);
hold off
legend('off')
grid on;
xlim([0 t(end)])

figure 
stem(1:iter, comp_time(1:iter) );
grid on
ylabel('Comp. time [s]')
xlabel('Time step')
xlim([0 iter])

%% Auxiliary function 
% Robot kinematics 
function [dq] = robot_kinematics(t, s, P, t0, tf, Omega_max)
    tau = 2 * (t-t0)/(tf-t0) - 1;
    tau = max(-1, min(1, tau));
    B = PolynomialBases.Legendre().basis(size(P,2)-1, tau);
    u = P * B;
    
%     idx2 = abs(u(1:4,:)) == max(abs(u(1:4,:)), [], 1);
%     idx = u(idx2) > Omega_max(1);
%     u(idx2(idx)) = repmat(Omega_max(1), sum(idx2(idx),1), 1);
%     
%     idx2 = abs(u(5:6,:)) == max(abs(u(5:6,:)), [], 1);
%     idx = u(idx2) > Omega_max(2);
%     u(idx2(idx)) = repmat(Omega_max(2), sum(idx2(idx),1), sum(idx));

    b = 1e-3; 
    k = 1e-3;
    dq = u + b * sin(t * u / 1E3) + k * (s(1:6,:) - s(1:6,:).^2);
end