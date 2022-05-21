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

function [C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_coe, final_coe, K, T, m, sampling_distribution, basis, n, setup)
    % Characteristics of the system 
    mu = system.mu;             % Characteristic gravitational parameter
    r0 = system.distance;       % Characteristic distance
    t0 = system.time;           % Characteristic time

    % Approximation order 
    n = repmat(n, [1 3]); 

    % Boundary conditions 
    % Initial state vector 
    s = coe2state(mu, initial_coe);
    initial = cylindrical2cartesian(s, false).';
    
    % Final state vector 
    s = coe2state(mu, final_coe);
    final = cylindrical2cartesian(s, false).';
    
    % Initial TOF
    tfapp = initial_tof(mu, T, initial, final);

    % Normalization
    % Gravitational parameter of the body
    mu = mu*(t0^2/r0^3);
    
    % Boundary conditions
    initial_coe(1) = initial_coe(1)/r0;
    final_coe(1) = final_coe(1)/r0;
    
    s = coe2state(mu, initial_coe);
    initial = cylindrical2cartesian(s, false).';
    
    s = coe2state(mu, final_coe);
    final = cylindrical2cartesian(s, false).';
    
    % Add additional revolutions 
    final(2) = final(2)+2*pi*K;
    
    % Time of flight
    tfapp = tfapp/t0;
    
    % Spacecraft propulsion parameters 
    T = T*(t0^2/r0);

    % Core optimization
    % Initial guess for the boundary control points
    mapp = 300;   
    tau = collocation_grid(mapp, sampling_distribution, '');
    [~, Capp, Napp, tfapp] = initial_approximation(tau, tfapp, initial, final, basis); 
    
    % Initial fitting for n+1 control points
    [P0, ~] = initial_fitting(n, tau, Capp, basis);
    
    % Final collocation grid and basis 
    tau = collocation_grid(m, sampling_distribution, '');
    [B, tau] = state_basis(n, tau, basis);

    % Initial guess 
    x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
    x0 = [x0; tfapp; Napp];
    L = length(x0)-2;
    
    % Upper and lower bounds (empty in this case)
    P_lb = [-Inf*ones(L,1); 0; 0];
    P_ub = [Inf*ones(L,1); Inf; Inf];
    
    % Objective function
    objective = @(x)cost_function(mu, initial, final, n, tau, x, B, basis, sampling_distribution);
    
    % Linear constraints
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    % Non-linear constraints
    nonlcon = @(x)constraints(mu, T, initial, final, n, x, B, basis, sampling_distribution, tau);
    
    % Modification of fmincon optimisation options and parameters (according to the details in the paper)
    options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'off', 'Algorithm', 'sqp');
    options.MaxFunctionEvaluations = 1e6;
    
    % Optimisation
    [sol, dV, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);
    
    % Solution 
    P = reshape(sol(1:end-2), [size(P0,1) size(P0,2)]);     % Optimal control points
    tf = sol(end-1);                                        % Optimal time of flight
    N = floor(sol(end));                                    % Optimal number of revolutions 
    
    P = boundary_conditions(tf, n, initial, final, N, P, B, basis);
    
    % Final state evolution
    C = evaluate_state(P,B,n);
    r = sqrt(C(1,:).^2+C(3,:).^2);
    
    % Solution normalization
    switch (sampling_distribution)
        case 'Sundman'
            % Normalised time grid
            options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);
            [s,time] = ode45(r,tau,0,options);
            s = tf*s; 
    
            % Control input
            u = acceleration_control(mu,C,tf,sampling_distribution);
            u = u/tf^2;
    
            % Trajectory cost
            dV = dV/tf;
    
        otherwise
            % Time domain normalization 
            switch (sampling_distribution)
                case 'Chebyshev'
                    tau = (1/2)*(1+tau);
                case 'Legendre'
                    tau = (1/2)*(1+tau);
                case 'Laguerre'
                    tau = collocation_grid(m, 'Legendre', '');
                    tau = (1/2)*(1+tau);
            end

            % Control input
            u = acceleration_control(mu,C,tf,sampling_distribution);
            u = u/tf^2;
    
            % Trajectory cost
            dV = dV/tf;
    end

    % Results 
    if (setup.resultsFlag)
        display_results(exitflag, output, r0, t0, tfapp, tf, dV);
        plots(system, tf, tau, C, u, T, initial_coe, final_coe, setup);
    end
end
 
