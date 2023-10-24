%% Project: SBOPT %%
% Date: 24/10/23

%% 3-DoF optimization of the chaser motion %% 
% This script provides a main interface to solve the general 3-DoF rotational problem for the robot %

% Inputs: - n, the vector of polynomial degrees
%         - m, the number of independent points in the grid
%         - S0, vector of initial conditions of the linear state 

% Output: - C, the optimal trajectory 
%         - P, the polynomial coefficients 
%         - u, the control trajectory
%         - omega, the angular velocity field

function [C, P, u, omega] = angular_command(n, m, S0)
    % Numerical solver definition 
    basis = 'Legendre';                   % Polynomial basis to be use
    time_distribution = 'Legendre';       % Distribution of time intervals
 
    solver = Solver(basis, n, time_distribution, m);

    % Problem definition
    Lc = 1;                               % Characteristic length [m]
    Tc = 300;                             % Characteristic time [s]
    Tmax = 0.004;                         % Maximum available torque [Nm]
    omega_max = 0.1;                      % Maximum angular velocity [rad/s]

    % Add attitude boundary conditions
    qf = [0.5;0.5;0.5;0.5];               % Final relative quaternion (null)
    omega_f = zeros(4,1);                 % Final relative angular velocity [rad/s]
    SF = [qf; omega_f];                   % Final conditions

    % Create the problem
    params(1) = Tc;                       % TOF [s]
    params(2) = Lc;                       % Maximum length [m]
    params(3) = Tmax;                     % Maximum control authority [Nm]
    params(4) = omega_max;                % Maximum rotational speed [rad/s]

% params(4) = mu;                  % Gauss constant
% params(5) = Orbit_t(2);          % Target orbital eccentricity
% params(6) = h;                   % Angular momentum magnitude
% params(7) = nu_0;                % Initial true anomaly [rad]
% params(8) = nu_f;                % Final true anomaly [rad]

    params(5:13) = diag([0.001 0.06 0.06]);    % Inertia tensor of the chaser [kg m^2]

    L = 2;                         % Degree of the dynamics (maximum derivative order of the ODE system)
    StateDimension = 4;            % Dimension of the configuration vector. Note the difference with the state vector
    ControlDimension = 3;          % Dimension of the control vector
    
    OptProblem = Problems.RobotRotoBerthing(S0, SF, L, StateDimension, ControlDimension, params);

    % Optimization    
    tic
    [C, ~, u, ~, ~, tau, ~, ~, P] = solver.solve(OptProblem);
    toc 

    % Angular velocity field
    omega = 2 * [C(4,:).*C(5,:)+C(3,:).*C(6,:)-C(2,:).*C(7,:)-C(1,:).*C(8,:); ...
                -C(3,:).*C(5,:)+C(4,:).*C(6,:)+C(1,:).*C(7,:)-C(2,:).*C(8,:); ...
                 C(2,:).*C(5,:)-C(1,:).*C(6,:)+C(4,:).*C(7,:)-C(3,:).*C(8,:)];

    % Final trajectory
    C = [tau; C(1:6,:)];
end