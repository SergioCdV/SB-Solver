%% Project: SBOPT %%
% Date: 01/08/22

%% 6-DoF optimization of the chaser motion %% 
% This script provides a main interface to solve the general 6-DoF problem for the robot %

% Inputs: - n, the vector of polynomial degrees
%         - m, the number of independent points in the grid
%         - S0, vector of initial conditions of the linear state 

% Output: - P, the optimal trajectory 
%         - C,  the polynomial coefficients
%         - u, the control traje
%         - omega, the angular velocity field
%         - maxIter, maximum number of iterations in the optimization
%         - P0, initial guess for the coefficientes

function [P, C, u, omega] = joint_command(n, m, S0, maxIter, P0)
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
    
    Tc = 300;                       % Characteristic time [s]
    Omega_max = [pi; 2*pi];         % Maximum angular velocity [rad/s]

    % Problem definition 
    L = 1;                          % Degree of the dynamics (maximum derivative order of the ODE system)
    StateDimension = 6;             % Dimension of the configuration vector. Note the difference with the state vector
    ControlDimension = 6;           % Dimension of the control vector
    
    % Linear problem data
    params(1) = 0;                  % Initial time [s]
    params(2) = Tc;                 % Final time [s]
    params(3:4) = Omega_max;        % Maximum control authority 
    
    % DH parameters of the robot
    SF = [0 -3*pi/4 +pi/2 -3*pi/4 pi/2 pi/2].';
    
    s_ref = [0.38 -0.1306 0.408 zeros(1,3) 1 1e-3 * rand(1,6)].';
    
    % Reference trajectory polynomial
    epsilon = 1e-3^2;                                                 % Numerical tolerance for the Jacobian determinant
    params(5) = epsilon;
    params(6:6+size(s_ref,1)-1) = reshape(s_ref(:,end), 1, []);  

    OptProblem = Problems.RobotDiffKinematics(S0, SF, L, StateDimension, ControlDimension, params);

    % Optimization    
    [P, ~, u, ~, ~, ~, ~, ~, C] = solver.solve(OptProblem);
end