%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 18/08/22
% File: mfls_optimization.m 
% Issue: 0 
% Validated: 18/08/22

%% Shapse-based optimization %%
% Function to compute the low-thrust orbital transfer in the CR3BP using a polynomial
% shape-based approach

% Inputs: - structure system, containing the physical information of the
%           CR3BP of interest
%         - vector initial_state, the initial parameters
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

function [C, dV, u, tf, tfapp, tau, exitflag, output] = mfls_optimization(wp, wv, kap, phi, psi, initial, T, setup)
    % Setup of the algorithm
    n = setup.order;                        % Order in the approximation of the state vector
    basis = setup.basis;                    % Polynomial basis to be used 
    sampling_distribution = setup.grid;     % Sampling grid to be used
    m = setup.nodes;                        % Number of nodes in the grid
    cost = setup.cost_function;             % Cost function to be minimized   

    % Approximation order 
    if (length(n) == 1)
        n = repmat(n, [1 3]);
    end

    % Initial TOF
    tfapp = 2*pi;

    % Final collocation grid and basis 
    tau = sampling_grid(m, sampling_distribution, '');
    [B, tau] = state_basis(n, tau, basis);

    % Boundary conditions   
    final = zeros(1,4);            % Rendezvous condition

    % Initial guess for the boundary control points
    mapp = 300;   
    tapp = sampling_grid(mapp, sampling_distribution, '');
    [~, Capp] = initial_approximation(tfapp, tapp, initial, final, basis); 
    
    % Initial fitting for n+1 control points
    [P0, ~] = initial_fitting(n, tapp, Capp, basis);
    
    % Initial guess reshaping
    x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
    L = length(x0);
    x0 = [x0; tfapp];
    
    % Upper and lower bounds 
    P_lb = [-Inf*ones(L,1); 0];
    P_ub = [Inf*ones(L,1); 10*12*2*pi];
    
    % Objective function
    objective = @(x)cost_function(wp, wv, kap, phi, psi, cost, initial, final, n, tau, x, B, basis);
    
    % Linear constraints and inequalities
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    % Non-linear constraints
    nonlcon = @(x)constraints(wp, wv, kap, phi, psi, tau, T, initial, final, n, x, B, basis);
    
    % Modification of fmincon optimisation options and parameters (according to the details in the paper)
    options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'off', 'Algorithm', 'sqp');
    options.MaxFunctionEvaluations = 1e6;
    
    % Optimisation
    [sol, dV, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);
    
    % Solution 
    P = reshape(sol(1:end-1), [size(P0,1) size(P0,2)]);   % Optimal control points
    tf = sol(end);                                        % Optimal time of flight
    
    % Final control points imposing boundary conditions
    P = boundary_conditions(tf, n, initial, final, P, B, basis);
    
    % Final state evolution
    C = evaluate_state(P,B,n);

    % Control input
    u = acceleration_control(wp, wv, kap, phi, psi, C, tf, tau);

    % Dimensionalization 
    C(3:4,:) = C(3:4,:)/tf;
    C(5:6,:) = C(5:6,:)/tf^2;

    % Time domain normalization 
    switch (sampling_distribution)
        case 'Chebyshev'
            tau = (1/2)*(1+tau);
            tf = tf*2;
        case 'Legendre'
            tau = (1/2)*(1+tau);
            tf = tf*2;
        case 'Laguerre'
            tau = collocation_grid(m, 'Legendre', '');
            tau = (1/2)*(1+tau);
            tf = tf*2;
    end
end