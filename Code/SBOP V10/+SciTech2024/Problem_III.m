%% Project: SBOPT %%
% Date: 18/01/23

%% Pre-capture phase %% 
% This script provides the solving of Problem III in URJC rendezvous laboratory,
% SciTech2024

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 10;                                % Polynomial order in the state vector expansion
m = 30;                                % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

Tc = 120;                              % Characteristic time [s]
Omega_max = [pi; 2*pi];                % Maximum angular velocity [rad/s]

%% Problem definition 
L = 1;                                 % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 6;                    % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 6;                  % Dimension of the control vector

% Linear problem data
params(1) = 0;                         % Initial time [s]
params(2) = 0.1 * Tc;                  % Final time [s]
params(3:4) = Omega_max;               % Maximum control authority 
params(5) = 5;                         % Constrain on the linear velocity
params(6) = 2;                         % Constraint on the angular velocity

epsilon = 1e-4^2;                      % Numerical tolerance for the Jacobian determinant
params(7) = epsilon;

% Final conditions
params(8:20) = [0.2; -0.05; 0.1; sin(deg2rad(0)/2) * [0;0;1]; cos(deg2rad(0)/2); zeros(3,1); deg2rad(0.1); deg2rad(-1); 0].';
params(8:20) = [0.4; -0.138; 0.38; sin(deg2rad(0)/2) * [0;0;1]; cos(deg2rad(0)/2); zeros(3,1); deg2rad(0.1); deg2rad(-0.5); 0].';

% DH parameters of the robot
S0 = [0 deg2rad(-151.31) deg2rad(142.22) deg2rad(-145) deg2rad(89.37) deg2rad(125)].';
%S0 = [0 deg2rad(-135) pi/2 deg2rad(-135) pi/2 pi/2].';
SF = [0 -3*pi/4 +pi/2 -3*pi/4 pi/2 0].';

%% Optimization (NMPC-RTI)
% Setup
options = odeset('AbsTol', 1e-22, 'RelTol', 2.25e-14);  % Integration tolerances
Ts = 1;                                                 % Sampling time [s]

% Noise 
Sigma_t = deg2rad(1)^2 * eye(6);                        % Joint state covariance [rad]

GoOn = true;                                            % Convergence boolean
iter = 1;                                               % Initial iteration
maxIter = 1e6;                                          % Maximum number of iterations
elapsed_time = 0;                                       % Elapsed time
tof = Tc;                                               % Maximum phase time

C = zeros(maxIter,12);
C(1,:) = [S0.' zeros(1,6)];
U = []; 
y0 = S0.';

while (GoOn && iter < maxIter)
    % Optimization (feedback phase)
    tic
    OptProblem = Problems.RobotDiffKinematics(S0, SF, L, StateDimension, ControlDimension, params);
    [T, dV, u, t0, tf, tau, exitflag, ~, P] = solver.solve(OptProblem);

    if (iter == 1)
        params(21:26) = mod(T(1:6,end), 2*pi);
        SF = mod(T(1:6,end), 2*pi);
        open_cost = dV;
        open_time = toc;
    else
        iter_time(iter-1) = toc;
    end
    
    u(1:3, max(abs(u(1:4,:))) > Omega_max(2)) = Omega_max(2);
    u(4:6, max(abs(u(4:6,:))) > Omega_max(2)) = Omega_max(2);

    % Preparing phase
    solver.InitialGuessFlag = true; 
    solver.P0 = P;
    solver.maxIter = 1;
    
    % Plant dynamics
    if (iter > 1)
        [~, s] = ode45(@(t,s)robot_kinematics(t, s, Pu, t0, tf), [Ts Ts + iter_time(iter-1)], y0, options);
        y0 = s(end,:);
    end

    Pu = PolynomialBases.Legendre().modal_projection(u);
    U = [U Pu * PolynomialBases.Legendre().basis(size(Pu,2)-1, 2* (linspace(0,Ts,100)-t0)/(tf-t0) -1)];     % Control trajectory

    [~, s] = ode45(@(t,s)robot_kinematics(t, s, Pu, t0, tf), [0 Ts], y0, options);  
    elapsed_time = (iter-1) * Ts;

    % Update initial conditions and state vector
    y0 = s(end,:);
    C(iter+1,:) = [s(end,:) u(:,1).'];

    % Navigation system
    S0 = C(iter+1,1:6).';
    S0 = S0 + mvnrnd(zeros(1,6), Sigma_t, 1).';    % Noisy state vector

    % Convergence 
    if (elapsed_time >= Tc)
        GoOn = false;
    else
        S0 = S0(1:6);

        % Update the number of iterations
        iter = iter + 1;
    end
end

U = [U zeros(size(U,1),1)];     % Control trajectory
t = Ts * (0:iter);              % Elapsed time vector
C = C(1:iter+1,:).';            % State trajectory

iter_time = mean(iter_time);

%% Analysis
% Check for singularities 
detJ = zeros(1,length(t)); 
v = zeros(6,length(t));

for i = 1:length(t)
    % Compute the Jacobian 
    [T, J] = Problems.RobotDiffKinematics.Kinematics(OptProblem.StateDim, ...
                                                     @(i,s)Problems.RobotDiffKinematics.ur3_dkinematics(OptProblem, i, s), ...
                                                     C(:,i));

    % Compute the determinant of the Jacobian
    detJ(i) = det(J);

    % Compute the frame velocities
    v(:,i) = J * C(7:12,i);
end

% Close loop cost 
close_cost = trapz(linspace(0, elapsed_time, size(U,2)), dot(U,U,1));

C(1:6,:) = mod(C(1:6,:), 2*pi);

%% Plots
% State representation
figure
hold on
xlabel('$t$ [s]')
ylabel('$\mathbf{\theta}$ [deg]')
plot(t, rad2deg(C(1:6,:)));
for i = 1:6
    yline( rad2deg(mod(params(20 + i), 2*pi)), 'k--')
end
legend('$\theta_1$', '$\theta_2$', '$\theta_3$', '$\theta_4$', '$\theta_5$', '$\theta_6$')
hold off
grid on;
xlim([0 t(end)])
yticklabels(strrep(yticklabels, '-', '$-$'));

% Propulsive acceleration plot
figure;
hold on
plot(linspace(0, elapsed_time, size(U,2)), max(abs(U(1:3,:)), [], 1), 'r');
plot(linspace(0, elapsed_time, size(U,2)), max(abs(U(4:6,:)), [], 1), 'b');
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
%% Auxiliary function 
% Robot kinematics 
function [dq] = robot_kinematics(t, s, P, t0, tf)
    tau = 2 * (t-t0)/(tf-t0) - 1;
    B = PolynomialBases.Legendre().basis(size(P,2)-1, tau);
    u = P * B;

    b = 1e-2; 
    k = 1e-3;
    dq = (1-b) * (u + b * sin(u.^2)) + k * (s(1:6,1) - s(1:6,1).^2);
end