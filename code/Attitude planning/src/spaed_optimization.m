%% Project: 
% Date: 19/05/22

%% Shapse-based optimization %%
% Function to compute the low-thrust orbital transfer using a polynomial
% shape-based approach

% Inputs: - structure system, containing the physical information of the
%           2BP of interest
%         - vector initial_bcs, the initial Tait-Bryan angles
%         - vector final_bcs, the final Tait-Bryan angles
%         - scalar T, the maximum allowed acceleration
%         - scalar m, the number of sampling nodes to use 
%         - string sampling_distribution, to select the sampling distribution
%           to use 
%         - string basis, the polynomial basis to be used in the
%           optimization
%         - scalar n, the polynomial degree to be used 
%         - structure setup, containing the setup of the figures

% Outputs: - array C, the final state evolution matrix
%          - scalar dV, the final dV cost of the transfer 
%          - array u, a 3xm matrix with the control input evolution 
%          - vector tau, the time sampling points final distribution
%          - exitflag, the output state of the optimization process 
%          - structure output, containing information on the final state of
%            the optimization process

function [C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_bcs, final_bcs, tf, T, m, sampling_distribution, basis, n, setup)
    % Characteristics of the system 
    I = system.Inertia;      % Inertia of the system

    % Approximation order 
    n = repmat(n, [1 4]); 

    % Boundary conditions 
    % Initial state vector 
    initial = TaitBryan2quat(false, initial_bcs);
    
    % Final state vector 
    final = TaitBryan2quat(false, final_bcs);
    
    % Initial TOF
    tfapp = tf;

    % Mapping to quaternion derivatives
    initial(5:8) = (1/2)*quaternion_product(initial(1:4), [0 initial(5:7)]);
    final(5:8) = (1/2)*quaternion_product(final(1:4), [0 final(5:7)]); 
    
    % Core optimization
    % Initial guess for the boundary control points
    mapp = 300;   
    tapp = collocation_grid(mapp, sampling_distribution, '');
    [~, Capp] = initial_approximation(tapp, tfapp, initial, final, basis); 
    
    % Initial fitting for n+1 control points
    [P0, ~] = initial_fitting(n, tapp, Capp, basis);
    
    % Final collocation grid and basis 
    tau = collocation_grid(m, sampling_distribution, '');
    [B, tau] = state_basis(n, tau, basis);

    % Initial guess 
    x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
    L = length(x0);
    
    % Upper and lower bounds (empty in this case)
    P_lb = [-Inf*ones(L,1)];
    P_ub = [Inf*ones(L,1)];
    
    % Objective function
    objective = @(x)cost_function(I, tf, initial, final, n, tau, x, B, basis);
    
    % Linear constraints
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    % Non-linear constraints
    nonlcon = @(x)constraints(system, tf, T, initial, final, n, x, B, basis);
    
    % Modification of fmincon optimisation options and parameters (according to the details in the paper)
    options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'off', 'Algorithm', 'sqp');
    options.MaxFunctionEvaluations = 1e6;
    
    % Optimisation
    [sol, dV, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);
    
    % Solution 
    P = reshape(sol, [size(P0,1) size(P0,2)]);                       % Optimal control points
    P = boundary_conditions(tf, n, initial, final, P, B, basis);
    
    % Final state evolution
    C = evaluate_state(P,B,n);
    
    % Solution normalization
    % Control input
    u = acceleration_control(I,C,tf);

    % Time domain normalization 
    switch (sampling_distribution)
        case 'Chebyshev'
            tau = (1/2)*(1+tau);
        case 'Legendre'
            tau = (1/2)*(1+tau);
            tf = tf*2;
        case 'Laguerre'
            tau = collocation_grid(m, 'Legendre', '');
            tau = (1/2)*(1+tau);
    end

    % Results 
    if (setup.resultsFlag)
        display_results(exitflag, output, tfapp, tf, dV);
        plots(system, tf, tau, C, u, T, setup);
    end
end