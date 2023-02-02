%% Project: 
% Date: 19/05/22

%% Shapse-based optimization %%
% Function to compute optimal trajectories

% Inputs: - vector initial, the initial state 
%         - vector final, the final state
%         - scalar tf, the time of flight
%         - structure setup, containing the setup of the algorithm in general

% Outputs: - array C, the final state evolution matrix
%          - scalar cost, the minimum cost function value
%          - array u, a 3xm matrix with the control input evolution  
%          - scalar tf, the final time of flight 
%          - vector tau, the time sampling points final distribution
%          - exitflag, the output state of the optimization process 
%          - structure output, containing information on the final state of
%            the optimization process

function [C, cost, u, tf, tau, exitflag, output] = sb_solver(initial, final, tf, setup)
    % Setup of the algorithm
    n = setup.order;                        % Order in the approximation of the state vector
    basis = setup.basis;                    % Polynomial basis to be used 
    sampling_distribution = setup.grid;     % Sampling grid to be used
    m = setup.nodes;                        % Number of nodes in the grid

    % Approximation order 
    if (length(n) == 1)
        n = repmat(n, [1 1]);
    end
         
    % Quadrature definition
    [tau, W, J] = quadrature(n, m, sampling_distribution);

    % Final state basis
    [B, tau] = state_basis(n, tau, basis);

    % Initial fitting for n+1 control points
    P0 = zeros(length(n), max(n)+1);
    P0 = boundary_conditions(tf, n, initial, final, P0, B, basis);

    % Initial guess reshaping
    x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
    M = length(x0);
    
    % Upper and lower bounds (ask on this!)
    P_ub = Inf*ones(M,1);
    P_lb = -Inf*ones(M,1); 

    % Objective function
    objective = @(x)cost_function(tf, initial, final, B, basis, n, tau, W, x);

    % Non-linear constraints
    nonlcon = @(x)constraints(tf, initial, final, B, basis, n, tau, x);
    
    % Linear constraints and inequalities
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    % Modification of fmincon optimisation options and parameters (according to the details in the paper)
    options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'off', 'Algorithm', 'sqp');
    options.MaxFunctionEvaluations = 1e6;
    
    % Optimisation
    [sol, cost, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);
    
    % Solution 
    P = reshape(sol, [size(P0,1) size(P0,2)]);     % Optimal control points
    
    % Final control points imposing boundary conditions
    P = boundary_conditions(tf, n, initial, final, P, B, basis);
    
    % Final state evolution
    C = evaluate_state(P,B,n);

    % Dimensional control input
    u = acceleration_control(C, tf, tau) / (tf*J)^2;

    % Dimensional velocity 
    C(2,:) = C(2,:)/(tf*J);
    C(3,:) = C(3,:)/(tf*J)^2;
    
    % Time domain normalization and scale preserving
    switch (sampling_distribution)
        case 'Chebyshev'
            tau = (1/2)*(1+tau);
        case 'Legendre'
            tau = (1/2)*(1+tau);
        otherwise
    end

    % Results 
    if (setup.resultsFlag)
        display_results(exitflag, output);
    end
end

