%% Project: SBOPT %%
% Date: 01/08/22

%% 6-DoF optimization of the chaser motion %% 
% This script provides a main interface to solve the general 6-DoF problem for the robot %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Bernstein';                    % Polynomial basis to be use
time_distribution = 'Bernstein';        % Distribution of time intervals
n = 7;                                 % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

Lc = 1;                         % Characteristic length [m]
Tc = 300;                  % Characteristic time [s]
Tmax = 0.004;                      % Maximum available torque [Nm]

%% Problem definition 
L = 2;                           % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 6;              % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;            % Dimension of the control vector

% Initial boundary conditions (robot-supplied)
S0 = zeros(L * StateDimension, 1);    

% Final boundary conditions (free)
SF = zeros(L * StateDimension, 1);

% Compute the reference trajectory 
s_ref = zeros((L+1) * StateDimension, m+1);

%% Create the problem
% Linear problem data
params(1) = 0;                   % TOF 
params(2) = 10;                  % Maximum length
params(3) = Tmax;                % Maximum control authority 
params(4:9) = zeros(6,1);        % All joints are revolute

% Reference trajectory polynomial
params(4:4+size(s_ref,1)*size(s_ref,2)-1) = reshape(s_ref, 1, []);  

OptProblem = Problems.RobotDiffKinematics(S0, SF, L, StateDimension, ControlDimension, params);

%% Optimization
% Simple solution    
tic
[C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
toc 

% Check for singularities 
detJ = zeros(1,length(tau)); 

for i = 1:length(tau)
    % Compute the Jacobian 

    % Compute the determinant of the Jacobian
    detJ(i) = 1;
end

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
plot(tau, u(1:3,:), 'LineWidth', 0.3)
plot(tau, max(abs(u), [], 1), 'k');
yline(Tmax, 'k--')
xlabel('$t$')
ylabel('$\mathbf{\tau}$')
legend('$\tau_x$', '$\tau_y$', '$\tau_z$', '$\|\mathbf{\tau}\|_\infty$', '$\tau_{max}$');
grid on;
xlim([0 tau(end)])

%% Matlab CSV
csvwrite("test_1.csv", [C; tau])