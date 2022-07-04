%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/04/22

%% MPC wrapper %%
% Function to compute the flow and control law needed to implement an MPC
% regulator based on shape-based optimization

% Inputs: - structure system, containing the physical information of the
%           2BP of interest
%         - vector initial_coe, the initial orbital elements 
%         - vector final_coe, the final orbital elements  
%         - scalar T, the maximum allowed acceleration
%         - string method, the optimal control policy to be used

% Outputs: - vector t, the time evolution
%          - array S, the final state evolution matrix
%          - array U, a 3xm matrix with the control input evolution  
%          - scalar dV, the final dV cost of the transfer 

function [t, S, U, dV] = MPC_wrapper(system, initial_coe, final_coe, T, method)            
    % Integration setup 
    RelTol = 2.25e-14; 
    AbsTol = 1e-22; 
    options = odeset('RelTol', RelTol, 'AbsTol', AbsTol);

    % Constants of the problem 
    mu = system.mu;                         % Gravitational parameter of the system
    system.mu = 1*mu;                       % True gravitational parameter of the system
    t0 = system.time;                       % Fundamental time unit
    r0 = system.distance;                   % Fundamental distance unit

    % Perturbations
    K = 0*[1e4/r0 1e4/r0 1e4/r0 1e1*(t0/r0) 1e1*(t0/r0) 1e1*(t0/r0)];
    
    % Initial integration    
    dt = 1e-2;                              % Time step (non-dimensional)

    % Initial conditions
    S = coe2state(1, initial_coe./[r0 ones(1,length(initial_coe)-1)]).';      
    
    % Setup of the MPC loop 
    GoOn = true;                            % Convergence boolean
    maxIter = 5e2;                          % Maximum number of division
    i = 1;                                  % Initial slice index

    % Setup of the shape-based method
    sampling_distribution = 'Regularized';  % Distribution of time intervals
    basis = 'Chebyshev';                    % Polynomial basis to be use
    n = 10;                                 % Polynomial order in the state vector expansion
    m = 100;                                % Number of sampling points
    setup.resultsFlag = false;              % Plotting options
    setup.animations = false;               % Animations options

    % Preallocation 
    t = [];                                  % Time of flight vector
    U = [];                                 % Control law vector
    
    % Main computation 
    while (GoOn && (i < maxIter))        
        % Compute the commands
        switch (method)
            case 'Shape-based'
                [~, ~, u, tf] = spaed_optimization(system, initial_coe, final_coe, 1, T, m, sampling_distribution, basis, n, setup);
            case 'Linear programming'
                [~, ~, u, tf] = spaed_optimization(system, initial_coe, final_coe, 0, T, m, sampling_distribution, basis, n, setup);
            case 'Non Linear programming'
                [~, ~, u, tf] = spaed_optimization(system, initial_coe, final_coe, 0, T, m, sampling_distribution, basis, n, setup);
            otherwise
                error('No valid control policy was selected.');
        end

        % Output saturation
        if (norm(u(:,1)) > T*(t0^2/r0))
            u(:,1) = zeros(3,1);
        else
            rho = cylindrical2cartesian(S(i,:).', false);
            R = [cos(rho(2)) sin(rho(2)) 0; -sin(rho(2)) cos(rho(2)) 0; 0 0 1].';
            u(:,1) = R*u(:,1);
        end

        % New integration time span (receeding horizon)
        tspan = 0:dt:tf;
                 
        % New integration
        [DT, St] = ode113(@(t,s)dynamics(1, 0, 0, t, s, u(:,1)), tspan, S(i,:), options);

        % True state evolution
        t = [t (i-1)*dt];
        S(i+1,:) = St(2,:);

        % Noisy evolution (simulating the navigation function)
        s = S(i,:) - K + 2.*K.*rand(1,6); 
        initial_coe = state2coe(mu/system.mu, s.', 'Inertial');
        initial_coe(1) = initial_coe(1)*r0;

        % Control law 
        U = [U u(:,1)];
        tf

        % Convergence 
        if (tf < dt)
            GoOn = false;
        else
            i = i+1; 
        end
    end

    % Normalization 
    dV = trapz(t,U)/tf*(r0/t0); % Final transfer cost
    t = t*t0;                   % Final time of flight
    U = U*(r0/t0^2);            % Final control law
end