%% Project: SBOPT %%
% Date: 01/08/22

%% 6-DoF optimization of the chaser motion %% 
% This script provides a main interface to solve the general 6-DoF problem for the robot %

% Inputs: - n, the vector of polynomial degrees
%         - m, the number of independent points in the grid
%         - S0, vector of initial conditions of the linear state 
%         - maxIter, maximum number of iterations in the optimization
%         - P0, initial guess for the coefficientes
%         - sf, the optimized boundary conditions

% Output: - P, the optimal trajectory 
%         - C,  the polynomial coefficients
%         - u, the control signal

function [P, C, u] = joint_command(n, m, S0, maxIter, P0, sf)
    % Numerical solver definition 
    basis = 'Legendre';                    % Polynomial basis to be use
    time_distribution = 'Legendre';        % Distribution of time intervals
 
    solver = Solver(basis, n, time_distribution, m);

    if (exist('P0', 'var'))
        solver.P0 = P0;
        solver.InitialGuessFlag = true;
    end

    if (exist('maxIter', 'var'))
        solver.maxIter = maxIter;
    end
    
    Tc = 120;                       % Characteristic time [s]
    Omega_max = [pi; 2*pi];         % Maximum angular velocity [rad/s]

    % Problem definition 
    L = 1;                          % Degree of the dynamics (maximum derivative order of the ODE system)
    StateDimension = 6;             % Dimension of the configuration vector. Note the difference with the state vector
    ControlDimension = 6;           % Dimension of the control vector
    
    % Linear problem data
    params(1) = 0;                  % Initial time [s]
    params(2) = 0.1 * Tc;           % Final time [s]
    params(3:4) = Omega_max;        % Maximum control authority 
    params(5) = 5;                  % Constrain on the linear velocity
    params(6) = 2;                  % Constraint on the angular velocity
    
    epsilon = 1e-4^2;               % Numerical tolerance for the Jacobian determinant
    params(7) = epsilon;

    % Final conditions
    params(8:20) = [0.2; -0.05; 0.1; sin(deg2rad(0)/2) * [0;0;1]; cos(deg2rad(0)/2); zeros(3,1); deg2rad(0.1); deg2rad(-1); 0].';
    
    if (exist('sf', 'var'))
        params(21:26) = mod(sf, 2*pi);
        SF = mod(sf, 2*pi);
    else
        SF = [0 -3*pi/4 +pi/2 -3*pi/4 pi/2 pi/2].';
    end

    % Create the problem
    OptProblem = Problems.RobotDiffKinematics(S0, SF, L, StateDimension, ControlDimension, params);

    % Optimization    
    [P, ~, u, ~, ~, ~, ~, ~, C] = solver.solve(OptProblem);
end