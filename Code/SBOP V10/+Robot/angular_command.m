%% Project: SBOPT %%
% Date: 24/10/23

%% 3-DoF optimization of the chaser motion %% 
% This script provides a main interface to solve the general 3-DoF rotational problem for the robot %

% Inputs: - n, the vector of polynomial degrees
%         - m, the number of independent points in the grid
%         - S0, vector of initial conditions of the linear state 

% Output: - P, the optimal trajectory 
%         - C,  the polynomial coefficients
%         - u, the control traje
%         - omega, the angular velocity field
%         - maxIter, maximum number of iterations in the optimization
%         - P0, initial guess for the coefficientes

function [P, C, u, omega] = angular_command(n, m, S0, maxIter, P0)
    % Numerical solver definition 
    basis = 'Legendre';                   % Polynomial basis to be use
    time_distribution = 'Legendre';       % Distribution of time intervals
 
    solver = Solver(basis, n, time_distribution, m);

    if (exist('P0', 'var'))
        solver.P0 = P0;
        solver.InitialGuessFlag = true;
    end

    if (exist('maxIter', 'var'))
        solver.maxIter = maxIter;
    end

    % Problem definition
    Lc = 1;                               % Characteristic length [m]
    Tc = 300;                             % Characteristic time [s]
    Tmax = 0.004;                         % Maximum available torque [Nm]
    omega_max = 0.1;                      % Maximum angular velocity [rad/s]

    % Add attitude boundary conditions
    S0 = S0.';
    q0 = S0(1:4);                                                      % Initial relative quaternion
    omega_0 = [S0(5:7); 0];                                            % Initial relative angular velocity of the chaser [rad/s]
    omega_0 = 0.5 * QuaternionAlgebra.right_isoclinic(omega_0) * q0;   % Quaternion kinematics
    S0 = [q0; omega_0];                                                % Initial conditions

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
    [P, ~, u, ~, ~, tau, ~, ~, C] = solver.solve(OptProblem);
    toc 

    % Angular velocity field
    omega = 2 * [P(4,:).*P(5,:)+P(3,:).*P(6,:)-P(2,:).*P(7,:)-P(1,:).*P(8,:); ...
                -P(3,:).*P(5,:)+P(4,:).*P(6,:)+P(1,:).*P(7,:)-P(2,:).*P(8,:); ...
                 P(2,:).*P(5,:)-P(1,:).*P(6,:)+P(4,:).*P(7,:)-P(3,:).*P(8,:)];

    % Final trajectory
    P = [tau; P(1:6,:)];
end