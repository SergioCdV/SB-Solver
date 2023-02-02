%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/04/22

%% MPC wrapper %%
% Function to compute the flow and control law needed to implement an MPC
% regulator based on shape-based optimization

% Inputs: - structure system, containing the physical information of the
%           rigig body of interest
%         - vector initial_state, the initial desired attitude stated  
%         - vector final_state, the final desired attitude state
%         - scalar T, the maximum allowed torque
%         - string method, the optimal control policy to be used
%         - structure setup, the general setup of the algorithm

% Outputs: - vector t, the time evolution
%          - array S, the final state evolution matrix
%          - array U, a 3xm matrix with the control input evolution  
%          - scalar dV, the final cost of the transfer 

function [t, S, U, dV] = MPC_wrapper(system, initial_state, final_state, T, method, setup)            
    % Integration setup 
    RelTol = 2.25e-14; 
    AbsTol = 1e-22; 
    options = odeset('RelTol', RelTol, 'AbsTol', AbsTol);

    % Characteristics of the system 
    I = system.Inertia;                     % Inertia matrix of the system

    % Perturbations
    K = [zeros(1,4) zeros(1,3)];
    
    % Initial integration    
    dt = 10;                                        % Time step (in second)

    % Initial conditions
    s = initial_state;                              % Observer state
    S = TaitBryan2quat(false, initial_state);       % Integration state

    qf = TaitBryan2quat(false, final_state);        % Final quaternion
    qf = quaternion_inverse(qf(1:4));               % Inverse of the final quaternion
    
    % Setup of the MPC loop 
    GoOn = true;                            % Convergence boolean
    maxIter = 5e2;                          % Maximum number of division
    i = 1;                                  % Initial slice index
    tol = 1e-3;                             % Convergence tolerance

    % Preallocation 
    t = zeros(1,maxIter);                   % Time of flight vector
    U = zeros(3,maxIter);                   % Control law vector

    % Initial guess for the TOF 
    tf = 1; 
    
    % Main computation 
    while (GoOn && (i < maxIter))        
        % Compute the commands
        switch (method)
            case 'Shape-based'
                [~, ~, u, tf] = spaed_optimization(system, s, final_state, 1, T, setup);
            case 'Linear programming'
                [~, ~, u, tf] = spaed_optimization(system, s, final_state, 1, T, setup);
            case 'Non Linear programming'
                [~, ~, u, tf] = spaed_optimization(system, s, final_state, 1, T, setup);
            otherwise
                error('No valid control policy was selected.');
        end

        % Output saturation
        if (norm(u(:,1)) > T)
            u(:,1) = zeros(3,1);
        end
                 
        % New integration
        [~, St] = ode113(@(t,s)attitude_dynamics(I, u(:,1), t, s), [0 dt], S(i,:), options);

        % True state evolution
        t(i) = (i-1)*dt;
        S(i+1,:) = St(end,:);

        % Noisy evolution (simulating the navigation function)
        s = S(i+1,:)+K.*normrnd(0,1,1,7); 
        s = TaitBryan2quat(true, s);

        % Control law 
        U(:,i) = u(:,1);

        % Convergence analysis
        error(1:4) = quaternion_product(S(i+1,1:4).', qf);                                   % Error to the final desired quaternion
        error(1) = 2*acos(error(1));                                                         % Error to the final desired quaternion
        error(2:4) = S(i+1,5:7)-final_state(4:6);                                            % Error to the final desired angular velocity
        if (norm(error) < tol)   
            GoOn = false;
        else
            i = i+1; 
        end
    end

    % Normalization 
    dV = trapz(t, sqrt(dot(U,U,1)));            % Final transfer cost
end

%% Auxiliary functions 
% Compute the derivative vectorfield of the attitude dynamics
function [ds] = attitude_dynamics(I,u,~,s)
    % State vector 
    q = s(1:4);         % Quaternion field 
    omega = s(5:7);     % Angular velocity field 

    % Euler equations 
    ds(1:4,1) = (1/2)*quaternion_product(q, [0; omega]);
    ds(5:7,1) = I^(-1)*(u-cross(omega, I*omega));
end