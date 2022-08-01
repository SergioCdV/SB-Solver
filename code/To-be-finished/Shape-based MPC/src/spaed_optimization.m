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
%         - scalar m, the number of sampling nodes to use 
%         - string dynamics, to indicate the parametrization of the
%           dynamics to be used
%         - string sampling_distribution, to select the sampling distribution
%           to use 
%         - string basis, the polynomial basis to be used in the
%           optimization
%         - scalar n, the polynomial degree to be used 
%         - structure setup, containing the setup of the figures

% Outputs: - array C, the final state evolution matrix
%          - scalar dV, the final dV cost of the transfer 
%          - array u, a 3xm matrix with the control input evolution  
%          - scalar tf, the final time of flight 
%          - scalar tfapp, the initial estimated time of flight 
%          - vector tau, the time sampling points final distribution
%          - exitflag, the output state of the optimization process 
%          - structure output, containing information on the final state of
%            the optimization process

function [C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_coe, final_coe, K, T, m, dynamics, sampling_distribution, basis, n, setup)
    % Characteristics of the system 
    mu = system.mu;             % Characteristic gravitational parameter
    r0 = system.distance;       % Characteristic distance
    t0 = system.time;           % Characteristic time

    % Approximation order 
    if (length(n) == 1)
        n = repmat(n, [1 3]);
    end

    % Boundary conditions 
    s = coe2state(mu, initial_coe);                     % Initial state vector 
    initial = cylindrical2cartesian(s, false).';        % Initial state vector in cylindrical coordinates
    s = coe2state(mu, final_coe);                       % Final state vector                   
    final = cylindrical2cartesian(s, false).';          % Final state vector in cylindrical coordinates 
    
    % Initial TOF
    tfapp = initial_tof(mu, T, initial, final);

    % Normalization
    mu = mu*(t0^2/r0^3);                                % Gravitational parameter of the body

    initial_coe(1) = initial_coe(1)/r0;                 % Boundary conditions normalization
    final_coe(1) = final_coe(1)/r0;                     % Boundary conditions normalization

    s = coe2state(mu, initial_coe);                     % Initial state vector 
    initial = cylindrical2cartesian(s, false).';        % Initial state vector in cylindrical coordinates
    s = coe2state(mu, final_coe);                       % Final state vector                   
    final = cylindrical2cartesian(s, false).';          % Final state vector in cylindrical coordinates 

    tfapp = tfapp/t0;                                   % Time of flight
    T = T*(t0^2/r0);                                    % Spacecraft propulsion parameters 

    % Add additional revolutions 
    final(2) = final(2)+2*pi*K;

    % Initial guess for the boundary control points
    mapp = 300;   
    tapp = sampling_grid(mapp, sampling_distribution, '');
    [~, Capp, Napp, tfapp] = initial_approximation(dynamics, tapp, tfapp, initial, final, basis); 
    
    % Initial fitting for n+1 control points
    [P0, ~] = initial_fitting(n, tapp, Capp, basis);
    
    % Final collocation grid and basis 
    tau = sampling_grid(m, sampling_distribution, 'Intersection');
    [B, tau] = state_basis(n, tau, basis);

    % Initial guess reshaping
    x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
    L = length(x0);
    x0 = [x0; tfapp; Napp];
    
    % Upper and lower bounds 
    P_lb = [-Inf*ones(L,1); 0; 0];
    P_ub = [Inf*ones(L,1); Inf; Inf];
    
    % Objective function
    objective = @(x)cost_function(mu, initial, final, n, tau, x, B, basis, dynamics);
    
    % Linear constraints and inequalities
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    % Non-linear constraints
    nonlcon = @(x)constraints(mu, T, initial, final, n, x, B, basis, dynamics);
    
    % Modification of fmincon optimisation options and parameters (according to the details in the paper)
    options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'off', 'Algorithm', 'sqp');
    options.MaxFunctionEvaluations = 1e6;
    
    % Optimisation
    [sol, dV, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);
    
    % Solution 
    P = reshape(sol(1:end-2), [size(P0,1) size(P0,2)]);     % Optimal control points
    tf = sol(end-1);                                        % Optimal time of flight
    N = floor(sol(end));                                    % Optimal number of revolutions 
    
    % Final control points imposing boundary conditions
    P = boundary_conditions(tf, n, initial, final, N, P, B, basis);
    
    % Final state evolution
    C = evaluate_state(P,B,n);

    % Control input
    u = acceleration_control(mu, C, tf, dynamics);
    u = u/tf^2;

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
    
    % Solution normalization
    switch (dynamics)
        case 'Regularized'
            % Initial TOF 
            rapp = sqrt(Capp(1,:).^2+Capp(3,:).^2);
            tfapp = tfapp*trapz(tapp, rapp);

            % Normalised time grid
            options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);
            [~, tau] = ode45(@(t,s)Sundman_transformation(basis, n, P, t, s), tau, 0, options);
    
            % Control input
            r = sqrt(C(1,:).^2+C(3,:).^2);
            u = u./(r.^2);
    
            % Final TOF 
            tf = tau(end)*tf;

        otherwise    
    end

    % Results 
    if (setup.resultsFlag)
        display_results(exitflag, output, r0, t0, tfapp, tf, dV);
        plots(system, tf, tau, C, u, T, initial_coe, final_coe, setup);
    end
end
 

%% Auxiliary functions 
% Compute the derivative of time with respect to the generalized anomaly 
function [dt] = Sundman_transformation(basis, n, P, t, s)
    B = state_basis(n,s,basis);
    C = evaluate_state(P,B,n);
    dt = sqrt(C(1,:).^2+C(3,:).^2);
end
