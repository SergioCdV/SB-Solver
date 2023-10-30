%% Project: SBOPT %%
% Date: 24/10/23

%% 3-DoF optimization of the chaser motion %% 
% This script provides a main interface to solve the general 3-DoF linear problem for the robot %

% Inputs: - n, the vector of polynomial degrees
%         - m, the number of independent points in the grid
%         - S0, vector of initial conditions of the linear state 

% Output: - P, the optimal trajectory 
%         - C,  the polynomial coefficients
%         - Cl, the polynomial coefficients in the physical space
%         - u, the control trajectory
%         - maxIter, maximum number of iterations in the optimization
%         - P0, initial guess for the coefficientes

function [P, C, Cl, u] = linear_command(n, m, S0, maxIter, P0)
    % Numerical solver definition 
    basis = 'Legendre';                                    % Polynomial basis to be use
    time_distribution = 'Legendre';                        % Distribution of time intervals

    solver = Solver(basis, n, time_distribution, m);

    if (exist('P0', 'var'))
        solver.P0 = P0;
        solver.InitialGuessFlag = true;
    end

    if (exist('maxIter', 'var'))
        solver.maxIter = maxIter;
    end
             
    % Problem definition
    mu = 3.986e14;                                         % Gravitational parameter of the Earth
    COE = [7011e3 0.00 deg2rad(190) deg2rad(98) 0 0];      % Target orbital elements
    h = sqrt(mu * COE(1) * (1-COE(2)^2));                  % Target specific angular momentum

    Lc = 1;                                                % Characteristic length [m]
    TOF = 100;                                             % Characteristic time [s]
    Fmax = 0.5e-1;                                         % Maximum available acceleration [m/s^2]

    n = sqrt(mu/COE(1)^3);                                 % Mean motion [rad/s]
    K = floor(TOF/(2*pi/n));                               % Number of complete revolutions
    dt = TOF - K * (2*pi/n);                               % Elapsed time in the last revolution [s]

    nu_0 = OrbitalDynamics.kepler(COE);                    % Initial true anomaly [rad]
    COE(6) = COE(6) + n * dt;                              % Final mean anomaly [rad]
    nu_f = 2*pi*K + OrbitalDynamics.kepler(COE);           % Final true anomaly [rad]
    
    % Final boundary conditions 
    omega = mu^2 / h^3;                                             % True anomaly angular velocity
    k = 1 + COE(2) * cos(nu_0);                                     % Transformation parameter
    kp =  - COE(2) * sin(nu_0);                                     % Derivative of the transformation
    L = [k * eye(3) zeros(3); kp * eye(3) eye(3)/(k * omega)];      % TH transformation matrix
    S0 = L * S0.';                                                  % TH initial boundary conditions
        
    % Assemble the state vector
    SF = [0.38; -0.1306; 0.408; zeros(3,1)];                        % Final conditions [m] [m/s]
    omega = mu^2 / h^3;                                             % True anomaly angular velocity
    k = 1 + COE(2) * cos(nu_f);                                     % Transformation parameter
    kp =  - COE(2) * sin(nu_f);                                     % Derivative of the transformation
    L = [k * eye(3) zeros(3); kp * eye(3) eye(3)/(k * omega)];      % TH transformation matrix
    SF = L * SF;                                                    % TH initial boundary conditions
    
    % Create the problem
    params(1) = TOF;                 % TOF [s]
    params(2) = Lc;                  % Maximum length
    params(3) = Fmax;                % Maximum control authority 
    params(4) = mu;                  % Gauss constant
    params(5) = COE(2);              % Target orbital eccentricity
    params(6) = h;                   % Angular momentum magnitude
    params(7) = nu_0;                % Initial true anomaly [rad]
    params(8) = nu_f;                % Final true anomaly [rad]
        
    L = 2;                           % Degree of the dynamics (maximum derivative order of the ODE system)
    StateDimension = 3;              % Dimension of the configuration vector. Note the difference with the state vector
    ControlDimension = 3;            % Dimension of the control vector
    
    OptProblem = Problems.RobotLinearBerthing(S0, SF, L, StateDimension, ControlDimension, params);
    
    % Optimization    
    tic
    [P, ~, u, ~, ~, tau, ~, ~, C] = solver.solve(OptProblem);
    toc
    
    % Dimensional space 
    for i = 1:length(tau) 
        omega = mu^2 / h^3;                                            % True anomaly angular velocity
        k = 1 + COE(2) * cos(tau(i));                                  % Transformation
        kp =  - COE(2) * sin(tau(i));                                  % Derivative of the transformation
        L = [k * eye(3) zeros(3); kp * eye(3) eye(3)/(k * omega)];     % TH transformation matrix
        P(1:6,i) = L \ P(1:6,i);                                       % Physical space
    end

     % Compute the elapsed time 
     cos_theta = cos(tau);
     sin_theta = sin(tau); 
     cos_E = (COE(2) + cos_theta) ./ (1 + COE(2) * cos_theta);
     sin_E = sqrt(1 - COE(2)^2) * sin_theta ./ (1 + COE(2) * cos_theta);
     E = atan2(sin_E, cos_E); 
     M = E - COE(2) * sin(E);
     tau = (M - M(1)) / n;
    
     % Final trajectory
     P = [tau; P(1:6,:)];

     % Compute the coefficients for the physical position trajectory
     [Cl] = PolynomialBases.Legendre().modal_projection(P(2:4,:));
     Cl = Cl(:,1:size(C,2));
end