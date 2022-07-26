%% Project: 
% Date: 19/05/22

%% Optical Orbit Determination %%
% Function to compute the orbital evolution of an optically observed object
% with an a-priori orbit and final orbit

% Inputs: - structure system, containing the physical information of the
%           2BP of interest
%         - scalar tf, the estimated epoch duration considered
%         - vector initial_coe, the initial orbital elements 
%         - vector final_coe, the final orbital elements 
%         - array measurements, an mx4 matrix of measurements in the form
%           of epoch | 3D vector measurement
%         - scalar T, the maximum allowed acceleration
%         - scalar m, the number of sampling nodes to use 
%         - scalar n, the polynomial degree to be used 
%         - string cost_policy, the strategy to minimize
%         - string dynamics, the parametrization of the dynamics
%           vectorfield to be used
%         - string sampling_distribution, to select the sampling distribution
%           to use 
%         - string basis, the polynomial basis to be used in the
%           optimization
%         - structure setup, containing the setup of the figures

% Outputs: - array C, the final state evolution matrix
%          - scalar e, the final e cost of the process
%          - array u, a 3xm matrix with the control input evolution  
%          - scalar tf, the final time of flight 
%          - scalar tfapp, the initial estimated time of flight 
%          - vector tau, the time sampling points final distribution
%          - exitflag, the output state of the optimization process 
%          - structure output, containing information on the final state of
%            the optimization process

function [C, e, u, tf, tfapp, tau, exitflag, output] = sbod_optimization(system, tf, initial_coe, final_coe, measurements, T, m, n, cost_policy, dynamics, sampling_distribution, basis, setup)
    % Sanity checks 
    if (isempty(T))
        T = Inf;
    end

    % Characteristics of the system 
    mu = system.mu;             % Characteristic gravitational parameter
    r0 = system.distance;       % Characteristic distance
    t0 = system.time;           % Characteristic time

    % Approximation order 
    if (length(n) == 1)
        n = repmat(n, [1 3]);
    end
    
    % Normalization
    mu = mu*(t0^2/r0^3);                            % Gravitational parameter of the body
    measurements(1,:) = measurements(1,:) / t0;     % Measurements epochs
    tfapp = tf/t0;                                  % Time of flight
    T = T*(t0^2/r0);                                % Spacecraft propulsion parameters 
    
    % Boundary conditions
    initial_coe(1) = initial_coe(1)/r0;    
    s = coe2state(mu, initial_coe);
    initial = cylindrical2cartesian(s, false).';
    
    if (~isempty(final_coe))
        % Final boundary conditions
        final_coe(1) = final_coe(1)/r0;
        s = coe2state(mu, final_coe);
        final = cylindrical2cartesian(s, false).';
    else
        final = [];
        basis = 'Bernstein'; 
        warning('No final boundary conditions have been specified')
    end

    % Initial guess for the boundary control points
    mapp = 300;   
    tapp =  sampling_grid(mapp, sampling_distribution, '');
    [~, Capp, Napp, tfapp] = initial_approximation(dynamics, tapp, tfapp, initial, final, basis); 
    
    % Initial fitting for n+1 control points
    [P0, ~] = initial_fitting(n, tapp, Capp, basis);
    
    % Final collocation grid and basis 
    tau = sampling_grid(m, sampling_distribution, '');
    [B, tau] = state_basis(n, tau, basis);

    % Initial guess 
    x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
    x0 = [x0; tfapp; Napp];
    L = length(x0)-2;
    
    % Upper and lower bounds (empty in this case)
    P_lb = [-Inf*ones(L,1); 0; 0];
    P_ub = [Inf*ones(L,1); Inf; Inf];
    
    % Objective function
    objective = @(x)cost_function(mu, initial, final, measurements, n, x, B, cost_policy, tau, sampling_distribution, basis, dynamics);
    
    % Linear constraints
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    % Non-linear constraints
    nonlcon = @(x)constraints(mu, T, initial, final, measurements, n, x, B, cost_policy, tau, sampling_distribution, basis, dynamics);
    
    % Modification of fmincon optimisation options and parameters (according to the details in the paper)
    options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'iter-detailed', 'Algorithm', 'sqp');
    options.MaxFunctionEvaluations = 1e6;
    
    % Optimisation
    [sol, e, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);
    
    % Solution 
    P = reshape(sol(1:end-2), [size(P0,1) size(P0,2)]);     % Optimal control points
    tf = sol(end-1);                                        % Optimal time of flight
    N = floor(sol(end));                                    % Optimal number of revolutions 
    
    P = boundary_conditions(tf, n, initial, final, N, P, B, basis);
    
    % Final state evolution
    C = evaluate_state(P,B,n);
    r = sqrt(C(1,:).^2+C(3,:).^2);

    % Control input
    u = acceleration_control(mu,C,tf,dynamics);
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
            [~, tau] = ode45(@(t,s)Sundman_transformation(basis,n,P,t,s), tau, 0, options);
    
            % Control input
            u = u./(r.^2);

            % Final cost function 
            switch (cost_policy)
                case 'Dynamic residual'
                    e = e / tf; 
            end
    
            % Final TOF 
            tf = tau(end)*tf;

        otherwise
            % Final cost function 
            switch (cost_policy)
                case 'Dynamic residual'
                    e = e / tf; 
            end
    end

    % Results 
    if (setup.resultsFlag)
        display_results(exitflag, output, r0, t0, tfapp, tf, cost_policy, e);
        plots(Capp, system, tf, tau, C, u, T, initial_coe, final_coe, setup);
    end
end
 

%% Auxiliary functions 
% Compute the derivative of time with respect to the generalized anomaly 
function [dt] = Sundman_transformation(basis, n, P, t, s)
    B = state_basis(n,s,basis);
    C = evaluate_state(P,B,n);
    dt = sqrt(C(1,:).^2+C(3,:).^2);
end
