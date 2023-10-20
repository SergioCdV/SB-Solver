%% Project: SBOPT %%
% Date: 01/08/22

%% 6-DoF optimization of the chaser motion %% 
% This script provides a main interface to solve the general 6-DoF problem for the robot %

%% Set up 
close all
clear

%% Compute some reference (optimal) linear and angular trajectories
% Linear evolution
Robot
s_ref = C(1:6,:);  

% Attitude evolution 
RobotRoto
s_ref = [s_ref(1:3,:); C(1:4,:); s_ref(4:6,:); omega];

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 7;                                % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

Lc = 1;                         % Characteristic length [m]
Tc = 300;                       % Characteristic time [s]
Omega_max = [pi; 2*pi];         % Maximum angular velocity [rad/s]

%% Problem definition 
L = 1;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 6;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 6;           % Dimension of the control vector

% Linear problem data
params(1) = 0;                  % Initial time [s]
params(2) = Tc;                 % Final time [s]
params(3:4) = Omega_max;        % Maximum control authority 

% DH parameters of the robot
S0 = [0 -pi/2 0 -pi/2 pi/2 0].';
SF = [0 -3*pi/4 +pi/2 -3*pi/4 pi/2 0].';

% Reference trajectory polynomial
params(5:5+size(s_ref,1)-1) = reshape(s_ref(:,end), 1, []);  

OptProblem = Problems.RobotDiffKinematics(S0, SF, L, StateDimension, ControlDimension, params);

%% Optimization
% Simple solution    
tic
[C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
toc 

% Check for singularities 
detJ = zeros(1,length(tau)); 

% for i = 1:length(tau)
%     % Compute the Jacobian 
%     [T] = Problems.RobotDiffKinematics.Kinematics(OptProblem.StateDim, ...
%                                                      @(i,s)Problems.RobotDiffKinematics.ur3_dkinematics(OptProblem, i, s), ...
%                                                      C(:,i));
% 
%     % Compute the determinant of the Jacobian
%     detJ(i) = det(J);
% end

% Dimensions
% C(1:3,:) = C(1:3,:) * Lc;
% tau = tau * 1000;

%% Plots
% State representation
figure
hold on
xlabel('$t$')
ylabel('$\mathbf{q}$')
plot(tau, C(1:6,:));
legend('$q_1$', '$q_2$', '$q_3$', '$q_4$', '$q_5$', '$q_6$')
hold off
legend('off')
grid on;
xlim([0 tau(end)])

figure
hold on
xlabel('$t$')
ylabel('$\mathbf{\omega}$')
plot(tau, C(7:12,:));
legend('$\omega_1$', '$\omega_2$', '$\omega_3$', '$\omega_4$', '$\omega_5$', '$\omega_6$')
hold off
legend('off')
grid on;
xlim([0 tau(end)])

figure
hold on
xlabel('$t$')
ylabel('det $J$')
plot(tau, detJ);
hold off
legend('off')
grid on;
xlim([0 tau(end)])

% Propulsive acceleration plot
figure;
hold on
plot(tau, u(1:6,:), 'LineWidth', 0.3)
plot(tau, max(abs(u(1:3,:)), [], 1), 'k');
plot(tau, max(abs(u(4:6,:)), [], 1), 'k');
yline(Omega_max(1), 'k--')
yline(Omega_max(2), 'k--')
xlabel('$t$')
ylabel('$\mathbf{\tau}$')
legend('$\tau_x$', '$\tau_y$', '$\tau_z$', '$\|\mathbf{\tau}\|_\infty$', '$\tau_{max,1}$', '$\tau_{max,2}$');
grid on;
xlim([0 tau(end)])

%% Matlab CSV
csvwrite("test_19.csv", [C; tau])