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
%         - vector uprev, the previously converged control profile
%         - structure setup, containing the setup of the algorithm in general

% Outputs: - array C, the final state evolution matrix
%          - scalar dV, the final dV cost of the transfer 
%          - array u, a 3xm matrix with the control input evolution  
%          - scalar tf, the final time of flight 
%          - scalar tfapp, the initial estimated time of flight 
%          - vector tau, the time sampling points final distribution
%          - scalar thetaf, the final true longitude
%          - exitflag, the output state of the optimization process 
%          - structure output, containing information on the final state of
%            the optimization process

function [C, dV, u, tf, tfapp, thetaf, tau, exitflag, output] = sb_solver(system, initial_coe, final_coe, TOF, T, setup)
    % Setup of the algorithm
    n = setup.order;                        % Order in the approximation of the state vector
    basis = setup.basis;                    % Polynomial basis to be used 
    sampling_distribution = setup.grid;     % Sampling grid to be used
    m = setup.nodes;                        % Number of nodes in the grid
    cost = setup.cost_function;             % Cost function to be minimized
    time_free = setup.FreeTime;             % Free or time-fixed problem

    % Characteristics of the system 
    mu = system.mu;                         % Characteristic gravitational parameter
    r0 = system.distance;                   % Characteristic distance
    t0 = system.time;                       % Characteristic time

    % Approximation order 
    if (length(n) == 1)
        n = repmat(n, [1 5]);
    end

    % Boundary conditions 
    s = coe2state(mu, initial_coe);                     % Initial state vector 
    initial = cylindrical2cartesian(s, false).';        % Initial state vector in cylindrical coordinates
    s = coe2state(mu, final_coe);                       % Final state vector                   
    final = cylindrical2cartesian(s, false).';          % Final state vector in cylindrical coordinates 
    
    % Initial TOF
    if (time_free)
        tfapp = initial_tof(mu, T, initial, final);
    else
        tfapp = max(TOF,0);
    end

    % Normalization
    gamma = r0/(t0/(1))^2;                              % Characteristic acceleration of the system
    mu = mu*(t0^2/r0^3);                                % Gravitational parameter of the body
    system.mu = mu; 

    initial_coe(1) = initial_coe(1)/r0;                 % Boundary conditions normalization
    initial = coe2equinoctial(initial_coe, true);       % Initial MEEs

    final_coe(1) = final_coe(1)/r0;                     % Boundary conditions normalization
    final = coe2equinoctial(final_coe, true);           % Final MEEs

    tfapp = tfapp/t0;                                   % Time of flight
    T = T/gamma;                                        % Spacecraft propulsion parameters 
 
    % Initial guess for the boundary control points
    mapp = 300;   
    tapp = sampling_grid(mapp, sampling_distribution, '');
    [~, Capp, thetaf] = initial_approximation(tapp, tfapp, initial, final, basis); 
    
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
    x0 = [x0; initial(end); thetaf; T];

    % Initial control guess 
    uprev = zeros(3,length(tau)); 
    
    % Upper and lower bounds 
    if (time_free)
        tol = 1e-8/gamma;
        P_lb = [-Inf*ones(L,1); -Inf; 0; T-tol];
        P_ub = [Inf*ones(L,1); Inf; Inf; T+tol];
    else
        P_lb = [-Inf*ones(L,1); -Inf; 0; 0];
        P_ub = [Inf*ones(L,1); Inf; Inf; 1/gamma];
    end

    % Longitude domain 
    switch (sampling_distribution)
        case 'Chebyshev'
            L = 0.5*(tau+1);
        case 'Legendre'
            L = 0.5*(tau+1);
        otherwise
            L = tau; 
    end
        
    % Linear constraints and inequalities
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    % Modification of fmincon optimisation options and parameters (according to the details in the paper)
    options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'off', 'Algorithm', 'sqp');
    options.MaxFunctionEvaluations = 1e6;
    
    % Loop 
    GoOn = true; 
    iter = 1; 
    maxIter = 2; 
    tol = 1e-4; 

    while (GoOn && iter < maxIter)
        % Objective function
        objective = @(x)cost_function(cost, mu, initial, final, uprev, B, basis, n, L, tau, W, x);
    
        % Non-linear constraints
        nonlcon = @(x)constraints(mu, initial, final, tfapp, time_free, uprev, B, basis, n, L, tau, x);
        
        % Optimisation
        [sol, dV, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);
        
        % Solution 
        P = reshape(sol(1:end-3), [size(P0,1) size(P0,2)]);     % Optimal control points
        theta0 = sol(end-2);                                    % Initial anomaly
        thetaf = sol(end-1);                                    % Final anomaly
        T = sol(end);                                           % Needed thrust vector
        
        % Final control points imposing boundary conditions
        P = boundary_conditions(n, initial(1:5), final(1:5), P, B, basis);
        
        % Final state evolution
        C = evaluate_state(P,B,n);
    
        % Compute the longitude evolution 
        L = theta0+thetaf*L;
    
        % Dimensional control input
        u = acceleration_control(mu, C, L, uprev) / thetaf;
     
        % Convergence analysis 
        ds = uprev-u;
    
        if (max(sqrt(dot(ds,ds,1))) < tol) 
            GoOn = false; 
        else
            % Next iteration variables 
            x0 = reshape(P, [size(P,1)*size(P,2) 1]);
            x0 = [x0; theta0; thetaf; T];
            uprev = u;
            iter = iter+1;
        end
    end

    % Dimensional velocity 
    C(6:10,:) = C(6:10,:) / thetaf;

    % Final time of flight 
    w = 1+C(2,:).*cos(L)+C(3,:).*sin(L);
    a = sqrt(mu*C(1,:)).*(w./C(1,:)).^2;
    for i = 1:size(C,2)
        B = control_input(mu, C(:,i)); 
        a(i) = 1/(a(i)+B(6,3)*u(3,i));
    end
    
    if (isempty(W))
        tf = thetaf*trapz(tau,a);
    elseif (length(W) ~= length(tau))
        tf = 0; 
        for i = 1:floor(length(tau)/length(W))
            tf = tf + thetaf*dot(W,a(1+length(W)*(i-1):length(W)*i));
        end
    else
        tf = thetaf*dot(W,a);
    end

    C = [C(1:5,:); L; C(6:10,:)];
    C = equinoctial2ECI(mu, C(1:6,:), true);

    % Domain normalization and scale preserving
    switch (sampling_distribution)
        case 'Chebyshev'
            tau = (1/2)*(1+tau);
            thetaf = thetaf/J;
            tfapp = tfapp/J;
        case 'Legendre'
            tau = (1/2)*(1+tau);
            thetaf = thetaf/J;
            tfapp = tfapp/J;
        otherwise
    end

    % Results 
    if (setup.resultsFlag)
        display_results(exitflag, cost, output, r0, t0, tfapp, tf, dV);
        plots(system, tf, tau, C, u, T, initial_coe, final_coe, setup);
    end
end

