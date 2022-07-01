%% Project: Shape-based attitude planning %%
% Date: 19/05/22

%% Shapse-based optimization %%
% Function to compute a constrained attitude path to satisfy given boundary
% conditions

% Inputs: - structure system, containing the physical information of the
%           rigid body of interest
%         - vector initial_bcs, the initial Tait-Bryan angles and angular
%           velocity
%         - vector final_bcs, the final Tait-Bryan angles and angular
%           velocity
%         - scalar T, the maximum allowed torque
%         - scalar m, the number of sampling nodes to use 
%         - string sampling_distribution, to select the sampling distribution
%           to use 
%         - string basis, the polynomial basis to be used in the
%           optimization
%         - scalar n, the polynomial degree to be used 
%         - structure setup, containing the setup of the figures

% Outputs: - array C, the final state evolution matrix
%          - scalar dV, the final domega cost of the transfer 
%          - array u, a 3xm matrix with the control input evolution 
%          - vector tau, the time sampling points final distribution
%          - exitflag, the output state of the optimization process 
%          - structure output, containing information on the final state of
%            the optimization process

function [C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_bcs, final_bcs, tf, T, m, sampling_distribution, basis, n, setup)
    % Characteristics of the system 
    I = system.Inertia;         % Inertia matrix of the system

    % Normalization constants 
    t0 = 3600;                 % Seconds in an hour

    % Approximation order 
    if (length(n) == 1)
        n = repmat(n, [1 4]);
    end

    % Boundary conditions 
    initial = TaitBryan2quat(false, initial_bcs);           % Initial state vector 
    final = TaitBryan2quat(false, final_bcs);               % Final state vector 

    % Normalization
    tfapp = tf/t0;                                          % Initial TOF
    T = T*t0^2;                                             % Maximum torque                                 
    
    % Mapping to quaternion derivatives
    initial(5:8) = (1/2)*quaternion_product(initial(1:4), [0 initial(5:7)]);
    final(5:8) = (1/2)*quaternion_product(final(1:4), [0 final(5:7)]); 
    
    % Initial guess for the boundary control points
    mapp = 300;   
    tapp = sampling_grid(mapp, sampling_distribution, '');
    [~, Capp] = initial_approximation(tapp, tfapp, initial, final, basis); 
    
    % Initial fitting for n+1 control points
    [P0, ~] = initial_fitting(n, tapp, Capp, basis);
    
    % Final collocation grid and basis 
    tau = sampling_grid(m, sampling_distribution, 'Intersection');
    [B, tau] = state_basis(n, tau, basis);

    % Initial guess 
    x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
    L = length(x0);
    x0 = [x0; tfapp];
    
    % Upper and lower bounds (empty in this case)
    P_lb = [-Inf*ones(L,1); 0];
    P_ub = [Inf*ones(L,1); 1];
    
    % Objective function
    objective = @(x)cost_function(I, initial, final, n, tau, x, B, basis);
    
    % Linear constraints
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    % Non-linear constraints
    nonlcon = @(x)constraints(system, T, initial, final, n, x, B, basis);
    
    % Modification of fmincon optimisation options and parameters (according to the details in the paper)
    options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'off', 'Algorithm', 'sqp');
    options.MaxFunctionEvaluations = 1e6;
    
    % Optimisation
    [sol, dV, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);
    
    % Solution 
    P = reshape(sol(1:end-1), [size(P0,1) size(P0,2)]);                   % Optimal control points
    tf = sol(end);                                                        % Final maneuver time

    P = boundary_conditions(tf, n, initial, final, P, B, basis);
    
    % Final state evolution
    C = evaluate_state(P,B,n);
    
    % Control input
    u = acceleration_control(I, C, tf);

    % Dimensional time units
    tf = tf*t0;
    tfapp = tfapp*t0;
    T = T/t0^2;
    u = u/tf^2;

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