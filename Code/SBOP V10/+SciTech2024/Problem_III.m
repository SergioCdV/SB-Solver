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
n = 7;                                 % Polynomial order in the state vector expansion
m = 50;                                % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

Tc = 120;                              % Characteristic time [s]
Omega_max = [pi; 2*pi];                % Maximum angular velocity [rad/s]

%% Problem definition 
L = 1;                                 % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 6;                    % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 6;                  % Dimension of the control vector

% Linear problem data
params(1) = 0;                         % Initial time [s]
params(2) = Tc;                        % Final time [s]
params(3:4) = Omega_max;               % Maximum control authority 

% DH parameters of the robot
S0 = [0 deg2rad(-135) pi/2 deg2rad(-135) pi/2 pi/2].';
SF = [0 -3*pi/4 +pi/2 -3*pi/4 pi/2 0].';

% Reference trajectory polynomial
epsilon = 1e-3^2;                           % Numerical tolerance for the Jacobian determinant
params(5) = epsilon;
params(6:11) = zeros(1,6);                  % Joints angular velocity  

%% Optimization (NMPC-RTI)
% Setup
options = odeset('AbsTol', 1e-22, 'RelTol', 2.25e-14);  % Integration tolerances
Ts = 30;                                                 % Sampling time

% Noise 
Sigma_t = 0.01 * eye(6);                                % Joint state covariance [rad]
Sigma_o = 0.05 * eye(6);                                % Joint velocity covariance [rad/s]

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
    elapsed_time = (iter-1) * Ts;
    OptProblem = Problems.RobotDiffKinematics(S0, SF, L, StateDimension, ControlDimension, params);
    [~, ~, u, t0, tf, tau, exitflag, ~, P] = solver.solve(OptProblem);
    
    if (norm(u(1:3,1), "inf") > Omega_max(1))
        u(norm(u(1:3,1), "inf") > Omega_max(1),1) = Omega_max(1);
    end
    if (norm(u(4:6,1), "inf") > Omega_max(2))
        u(norm(u(4:6,1), "inf") > Omega_max(2),1) = Omega_max(2);
    end
    U = [U u(:,1)];     % Control trajectory

    % Preparing phase
    solver.InitialGuessFlag = true; 
    solver.P0 = P;
    solver.maxIter = 1;
    
    % Plant dynamics 
    [~, s] = ode45(@(t,s)robot_kinematics(t, s, U(:,end)), linspace(0, Ts, 100), y0, options);  

    % Update initial conditions and state vector
    y0 = s(end,:);
    C(iter+1,:) = [s(end,:) u(:,1).'];

    % Navigation system
    S0 = C(iter+1,:);
    S0 = mvnrnd(S0, blkdiag(Sigma_t, Sigma_o), 1).';    % Noisy state vector

    % Convergence 
    if (elapsed_time >= Tc)
        GoOn = false;
    else
        % Shrinking horizon
        if (elapsed_time > 0.0 * Tc)
            tof = tof - Ts; 
            params(2) = tof;
        end

        params(6:11) = S0(7:12);                         % Joints angular velocity
        S0 = S0(1:6);

        % Update the number of iterations
        iter = iter + 1;
    end
end

U = [U zeros(size(U,1),1)];     % Control trajectory
t = Ts * (0:iter);              % Elapsed time vector
C = C(1:iter+1,:).';            % State trajectory

%% Analysis
% Check for singularities 
detJ = zeros(1,length(t)); 

for i = 1:length(t)
    % Compute the Jacobian 
    [T, J] = Problems.RobotDiffKinematics.Kinematics(OptProblem.StateDim, ...
                                                     @(i,s)Problems.RobotDiffKinematics.ur3_dkinematics(OptProblem, i, s), ...
                                                     C(i,:).');

    % Compute the determinant of the Jacobian
    detJ(i) = det(J);
end

%% Plots
% State representation
figure
hold on
xlabel('$t$ [s]')
ylabel('$\mathbf{\theta}$ [rad]')
plot(t, C(1:6,:));
legend('$\theta_1$', '$\theta_2$', '$\theta_3$', '$\theta_4$', '$\theta_5$', '$\theta_6$')
hold off
legend('off')
grid on;
xlim([0 t(end)])

figure
hold on
xlabel('$t$')
ylabel('$\dot{\mathbf{\theta}}$ [rad/s]')
plot(t, C(7:12,:));
legend('$\dot{\theta}_1$', '$\dot{\theta}_2$', '$\dot{\theta}_3$', '$\dot{\theta}_4$', '$\dot{\theta}_5$', '$\dot{\theta}_6$')
hold off
legend('off')
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

% Propulsive acceleration plot
figure;
hold on
plot(t, max(abs(U(1:3,:)), [], 1), 'r');
plot(t, max(abs(U(4:6,:)), [], 1), 'b');
yline(Omega_max(1), 'r--')
yline(Omega_max(2), 'b--')
xlabel('$t$')
ylabel('$\mathbf{\omega}$')
legend('$\|\omega^{1:3}\|_{\infty}$', '$\|\omega^{4:6}\|_{\infty}$', '$\tau_{max,1}$', '$\tau_{max,2}$');
grid on;
xlim([0 tau(end)])

%% Auxiliary function 
% Robot kinematics 
function [dq] = robot_kinematics(t,s,u)
    dq = u;
end