%% Project: SBOPT %%
% Date: 13/04/24

%% Attitude slew with RW %% 
% This script provides a main interface to solve an attitude slew problem usign Reaction Wheels %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 10;                        % Polynomial order in the state vector expansion
m = 500;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 2;                                % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 6;                   % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;                 % Dimension of the control vector

% Problem parameters
t0 = 0;                                % Initial clock 
tf = 120;                              % Final clock
I = [0.6 0.1 0; 0.1 1 0.2; 0 0.2 1] * 1E-2;    % Inertia dyadic
T = 0.0004;                             % Maximum torque [Nm]
omega_max = 0.1;                       % Maximum angular velocity [rad/s]
h_max = 0.002;                         % Maximum angular momentum [Nms]

problem_params = [t0 tf T omega_max h_max reshape(I, [], 9)];

% Add boundary conditions
S0 = [zeros(1,12)].';
SF = [1/3 1/3 1/3 zeros(1,9)].';

% Create the problem
OptProblem = Problems.AttitudeSlew(S0, SF, L, StateDimension, ControlDimension, problem_params);

%% Optimization
% Simple solution    
tic
[C, dV, u, t0, tf, tau, exitflag, output, P] = solver.solve(OptProblem);
toc 

P = PolynomialBases.Legendre().modal_projection(u(1:3,:));
u_con = @(t)Controller(L, n, basis, t0, tf, P, problem_params, t);

% Average results 
iter = 0; 
time = zeros(1,iter);
setup.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
    time(i) = toc;
end

%% Feasibility analysis 
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1E-22);
tspan = [t0, tf];
[t, s] = ode45(@(t,s)AttitudeDynamics(I, u_con, t, s), tspan, S0(1:9,1), options);

%% Plots
% State representation
figure_orbits = figure;
hold on
xlabel('$t$')
ylabel('$\mathbf{\sigma}$')
plot(tau, C(1:3,:));
plot(t, s(:,1:3), 'ok')
hold off
legend('$\sigma_x$', '$\sigma_y$', '$\sigma_z$')
grid on;
yticklabels(strrep(yticklabels, '-', '$-$'));

figure;
hold on
xlabel('$t$')
ylabel('$\mathbf{h}$')
plot(tau, C(4:6,:));
plot(tau, max(abs(C(4:6,:))))
yline(h_max, 'k--');
hold off
legend('$h_x$', '$h_y$', '$h_z$', '$\|\mathbf{h}\|_{\infty}$', '$h_{max}$')
grid on;
yticklabels(strrep(yticklabels, '-', '$-$'));

% Torque plot
figure_propulsion = figure;
hold on
plot(tau, u)
plot(tau, sqrt(dot(u,u,1)))
%yline(T, 'k--');
xlabel('$t$')
ylabel('$\mathbf{\tau}$')
legend('$\tau_x$', '$\tau_y$', '$\tau_z$', '$\|\mathbf{\tau}\|_2$', '$\tau_{max}$')
grid on;
yticklabels(strrep(yticklabels, '-', '$-$'));

%% Auxiliary function 
function [ds] = AttitudeDynamics(I, Controller, t, s)
    % Variables 
    sigma = s(1:3,:);       % MRPs
    omega = s(4:6,:);       % Angular velocity
    h = s(7:9,:);           % Angular momentum

    idx = dot(sigma, sigma, 1) > 1;
    sigma(:,idx) = -sigma(:,idx) ./ dot(sigma(:,idx), sigma(:,idx), 1);
    q = [sigma; -ones(1,size(sigma,2))];                                    % Modified MRPs

    % Kinematics
    ds(1:3,:) = 0.25 * QuaternionAlgebra.Quat2Matrix( q ) * omega; 

    % Compute the control law 
    u = Controller(t);

    % Dynamics 
    ds(4:6,:) = I \ (u + cross(I * omega, omega));  % Angular velocity dynamics
    ds(7:9,:) = -u + cross(h, omega);               % Angular momentum dynamics
end

% Control function 
function [u] = Controller(L, n, basis, t0, tf, P, params, t)
%     % Compute polynomial time 
    tau = 2 * (t - t0) / (tf-t0) - 1;
    tau = min(tau, 1);
    tau = max(tau, -1);
% 
%     % Evaluate the state function
%     B = Solver().state_basis(L, repmat(n, 1, 6), basis, tau);
%     s = Solver().evaluate_state(repmat(n, 1, 6), L, P, B);
% 
%     % Compute the associate control function 
%     % Constants
%     I = reshape(params(6:end), [3 3]);      % Inertia matrix
% 
%     % Pre-allocation
%     omega = zeros(3, size(tau,2));          % Angular velocity
%     alpha = omega;                          % Angular acceleration
%     sigma = s(1:3,:);
%     dsigma = s(7:9,:);
%     ddsigma = s(13:15,:);
% 
%     % Shadow transformation
%     idx = dot(sigma, sigma, 1) > 1;
%     sigma(:,idx) = -sigma(:,idx) ./ dot(sigma(:,idx), sigma(:,idx), 1);
% 
%     q = [sigma; -ones(1,size(tau,2))];      % Modified MRPs
%     q_squared = dot( q, q ).^2;             % Dot product of the associated quaternions
%   
%     omega = dsigma ./ q_squared;
%     
%     for i = 1:size(tau,2)
%         B = QuaternionAlgebra.Quat2Matrix( q(:,i) ).';
%         omega(:,i) = 4 * B * omega(:,i);
% 
%         dB = QuaternionAlgebra.dBmatrix( sigma(:,i), omega(:,i) );
%         alpha(:,i) = 4 * B * ( ddsigma(:,i) - 0.25 * dB * omega(:,i) );
%     end
% 
%     alpha = alpha ./ q_squared;                                                 % Angular acceleration
%     
%     % Torque solution
%     u = I * alpha + cross(omega, I * omega);

    B = PolynomialBases.Legendre().basis(size(P,2)-1, tau);
    u = P * B;
end