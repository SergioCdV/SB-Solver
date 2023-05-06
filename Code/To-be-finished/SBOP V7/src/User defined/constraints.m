%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - scalar mu, the gravitational parameter of the central body 
%         - scalar T, the maximum acceleration allowed for the spacecraft
%         - vector initial, the initial boundary conditions of the
%           trajectory 
%         - vector final, the initial boundary conditions of the
%           trajectory
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - vector x, the degree of freedom to be optimized 
%         - cell array B, the polynomial basis to be used
%         - string basis, the polynomial basis to be used
%         - string dynamics, the independent variable parametrization to be
%           used
%         - string cost, the cost function to be minimized

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(cost, mu, initial, final, B, basis, n, tau, x)
    % Extract the optimization variables
    P = reshape(x(1:end-9), [length(n), max(n)+1]);     % Control points
    lambda0 = x(end-8:end-3);                           % Initial co-state guesses
    tf = x(end-2);                                      % Final time of flight 
    thetaf = x(end-1);                                  % Final insertion phase
    T = x(end);                                         % Needed thrust vector

    % Boundary conditions points
    P = boundary_conditions(tf, n, initial, final, thetaf, P, B, basis);

    % Trajectory evolution
    C = evaluate_state(P,B,n);

    % Control input 
    [u, ~] = acceleration_control(mu, C, tf);

    % Equalities 
    ceq = [cos(thetaf)-cos(final(2)) sin(thetaf)-sin(final(2))];

    switch (cost)
        case 'Minimum time'
            % Integrate the co-states 
            lambda = repmat(lambda0.', size(C,2), 1);
            dynamics = @(tau, x)dual_dynamics(tf, C, x);
            [lambda, ~] = MCPI(tau, lambda, dynamics, length(tau)-1, 1e-5);
    
            % Compute the Hamiltonian constraint
            H = lambda(end,:)*C(4:9,end)+1;

        case 'Minimum fuel'    
            % Integrate the co-states 
            lambda = repmat(lambda0.', size(C,2), 1);
            dynamics = @(tau, x)dual_dynamics(tf, C, x);
            [lambda, ~] = MCPI(tau, lambda, dynamics, length(tau)-1, 1e-5);

            % Compute the primer vector constraint
            Lambda = lambda(:,4:6)./sqrt(dot(lambda(:,4:6),lambda(:,4:6),2));
            U = u./sqrt(dot(u,u,1));
            H = 1+dot(Lambda, U.', 1);

        otherwise
    end

    ceq = [ceq H];

    % Inequalities
    c = [u(1,:).^2+u(2,:).^2+u(3,:).^2-(tf^2*repmat(T,1,size(u,2))).^2];
end

%% Auxiliary functions 
% Dual dynamics 
function [dlambda] = dual_dynamics(tf, C, lambda)
    % Preallocation 
    dlambda = zeros(size(C,2),6); 

    A = [zeros(3) eye(3); zeros(3,6)];

    for i = 1:size(C,2)
        r = norm(C(1:3,i));
        A(4:6,1:3) = [-1/r^3+3*C(1,i)^2/r^5 3*C(1,i)*C(2,i)/r^5 3*C(1,i)*C(3,i)/r^5; ...
                      3*C(1,i)*C(2,i)/r^5 -1/r^3+3*C(2,i)^2/r^5 3*C(2,i)*C(3,i)/r^5; ...
                      3*C(1,i)*C(3,i)/r^5 3*C(2,i)*C(3,i)/r^5 -1/r^3+3*C(3,i)^2/r^5]; 
        dlambda(i,:) = -tf*lambda(i,:)*A;
    end
end

% Modified-Chebyshev Picard Iterations
function [x, state] = MCPI(tau, x0, dynamics, order, tol)
    % Set up of the method 
    GoOn = true;                            % Convergence boolean
    iter = 1;                               % Initial iteration
    maxIter = 1e2;                          % Maximum number of iterations
    error = 1e3;                            % Initial error

    % Constants 
    N = order+1;                            % Number of approximation terms

    T = CH_basis('first', order, tau).';

    W = diag(ones(1,N));                    % Weight matrix for the state vector
    W(1,1) = 1/2;                           % Weight matrix for the state vector

    V = (2/N)*diag(ones(1,N));              % Approximation weights for the dynamics
    V(1,1) = 1/N;                           % Approximation weights for the dynamics
    V(N,N) = 1/N;                           % Approximation weights for the dynamics

    R = (1/2)./(1:N-1);                     % Approximation weights for the dynamics
    R = diag([1 R]);                        % Approximation weights for the dynamics

    % Chebyshev polynomials difference
    S = zeros(N,N);
    S(1,1) = 1; 
    S(1,2) = -1/2;
    S(1,N) = (-1)^(order+1)/(order-1);
    for i = 2:order-1
        S(1,i+1) = (-1)^(i+1)*(1/(i-1)-1/(i+1));
        S(i,i-1) = 1; 
        S(i,i+1) = -1;
    end
    S(order,end-2) = 1;
    S(order,end) = -1;
    S(N,end-1) = 1;

    chi = zeros(N,size(x0,2));              % Boundary conditions term

    Ca = R*S*T.'*V;                         % Dynamics approximation matrix
    Cx = T*W;                               % State approximation matrix

    % Main loop
    while (GoOn && iter < maxIter)
        % Evaluate the dynamics 
        g = dynamics(tau,x0); 

        % Compute the beta coefficients
        chi(1,:) = 2*x0(1,:);
        beta = Ca*g+chi; 

        % Compute the new state trajectory
        x = Cx*beta;

        % Convergence check 
        dX = x-x0; 
        dX = sqrt(dot(dX,dX,2));
        if (error < tol && norm(dX) < tol)
            GoOn = false;
        else
            iter = iter+1;
            x0 = x; 
            error = norm(dX);
        end
    end

    % Final output 
    state.State = ~GoOn; 
    state.Iterations = iter; 
    state.Error = error; 
end