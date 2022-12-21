%% Project: 
% Date: 19/05/22

%% Shapse-based optimization %%
% Function to compute the low-thrust orbital transfer using a polynomial
% shape-based approach

% Inputs: - structure system, containing the physical information of the
%           2BP of interest
%         - vector initial_coe, the initial orbital elements 
%         - vector final_coe, the final orbital elements 
%         - scalar K, an initial desired revolutions value 
%         - scalar T, the maximum allowed acceleration 
%         - structure setup, containing the setup of the algorithm in general

% Outputs: - array C, the final state evolution matrix
%          - scalar dV, the final dV cost of the transfer 
%          - array u, a 3xm matrix with the control input evolution  
%          - scalar tf, the final time of flight 
%          - scalar tfapp, the initial estimated time of flight 
%          - vector tau, the time sampling points final distribution
%          - exitflag, the output state of the optimization process 
%          - structure output, containing information on the final state of
%            the optimization process

function [C, dV, u, tf, tfapp, tau, exitflag, output] = sb_solver(system, initial, final, TOF, T, setup)
    % Setup of the algorithm
    n = setup.order;                        % Order in the approximation of the state vector
    basis = setup.basis;                    % Polynomial basis to be used 
    sampling_distribution = setup.grid;     % Sampling grid to be used
    m = setup.nodes;                        % Number of nodes in the grid
    cost = setup.cost_function;             % Cost function to be minimized
    time_free = setup.FreeTime;             % Free or time-fixed problem

    % Characteristics of the system 
    M1 = system.M1;                         % Characteristic mass of the first primary
    M2 = system.M2;                         % Characteristic mass of the second primary
    r0 = system.distance;                   % Characteristic distance
    t0 = system.time;                       % Characteristic time

    % Approximation order 
    if (length(n) == 1)
        n = repmat(n, [1 3]);
    end
    
    % Normalization
    mu = M2/(M1+M2);                                    % Gravitational parameter of the body
    system.mu = mu; 
    
    initial = [initial(1:3)/r0 initial(4:6)*(t0/r0)];  % Initial state vector in Cartesian coordinates
    final = [final(1:3)/r0 final(4:6)*(t0/r0)];        % Final state vector in Cartesian coordinates

    gamma = r0/(t0/(1))^2;                              % Characteristic acceleration of the system
    T = T/gamma;                                        % Spacecraft propulsion parameters 

    % Initial TOF
    if (time_free)
        tfapp = initial_tof(mu, T, initial, final) / t0;
    else
        tfapp = max(TOF,0) / t0;
    end
    
    tfapp = 2*pi; 

    % Transformation to synodic coordinates 
    initial(4:6) = initial(4:6)+cross([0;0;1], initial(1:3));
    final(4:6) = final(4:6)+cross([0;0;1], final(1:3));
 
    % Initial guess for the boundary control points
    mapp = 300;   
    tapp = sampling_grid(mapp, sampling_distribution, '');
    [~, Capp, thetaf, tfapp] = initial_approximation(tapp, tfapp, initial, final, basis); 
    
    % Initial fitting for n+1 control points
    [P0, ~] = initial_fitting(n, tapp, Capp, basis);
    
    % Quadrature definition
    [tau, W, J] = quadrature(n, m, sampling_distribution);

    % Final TOF scaling
    tfapp = tfapp*J;

    % Final state basis
    [B, tau] = state_basis(n, tau, basis);

    % Initial guess reshaping
    x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
    L = length(x0);
    x0 = [x0; tfapp; thetaf; T];
    
    % Upper and lower bounds 
    if (time_free)
        tol = 1e-8/gamma;
        P_lb = [-Inf*ones(L,1); 0; 0; T-tol];
        P_ub = [Inf*ones(L,1); Inf; Inf; T+tol];
    else
        tol = 1e-4;
        P_lb = [-Inf*ones(L,1); max(tfapp-tol,0); 0; 0];
        P_ub = [Inf*ones(L,1); tfapp+tol; Inf; 1/gamma];
    end
    
    % Objective function
    objective = @(x)cost_function(cost, mu, initial, final, B, basis, n, tau, W, x);

    % Non-linear constraints
    nonlcon = @(x)constraints(mu, initial, final, B, basis, n, tau, x);
    
    % Linear constraints and inequalities
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    % Modification of fmincon optimisation options and parameters (according to the details in the paper)
    options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'off', 'Algorithm', 'sqp');
    options.MaxFunctionEvaluations = 1e6;
    
    % Optimisation
    [sol, dV, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);
    
    % Solution 
    P = reshape(sol(1:end-3), [size(P0,1) size(P0,2)]);     % Optimal control points
    tf = sol(end-2);                                        % Optimal time of flight
    thetaf = sol(end-1);                                    % Final anomaly
    T = sol(end);                                           % Needed thrust vector
    
    % Final control points imposing boundary conditions
    P = boundary_conditions(tf, n, initial, final, thetaf, P, B, basis);
    
    % Final state evolution
    C = evaluate_state(P,B,n);

    % Dimensional control input
    u = acceleration_control(mu, C, tf) / tf^2;

    % Dimensional velocity and acceleration
    C(4:6,:) = C(4:6,:)/tf;
    C(7:9,:) = C(7:9,:)/tf^2;
    
    % Time domain normalization and scale preserving
    switch (sampling_distribution)
        case 'Chebyshev'
            tau = (1/2)*(1+tau);
            tf = tf/J;
            tfapp = tfapp/J;
        case 'Legendre'
            tau = (1/2)*(1+tau);
            tf = tf/J;
            tfapp = tfapp/J;
        otherwise
    end

    % Results 
    if (setup.resultsFlag)
        display_results(exitflag, cost, output, r0, t0, tfapp, tf, dV);
        plots(system, tf, tau, C, u, T, initial, final, setup);
    end
end

