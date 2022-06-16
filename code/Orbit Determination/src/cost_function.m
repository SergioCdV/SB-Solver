%% Project: 
% Date: 31/01/22

%% Cost function %%
% Function to estimate the time of flight

% Inputs: - scalar mu, the gravitational parameter of the central body 
%         - vector initial, the initial boundary conditions of the
%           trajectory 
%         - vector final, the initial boundary conditions of the
%           trajectory
%         - array measurements, an mx4 matrix of measurements in the form
%           of epoch | 3D vector measurement
%         - vector x, the degree of freedom to be optimized 
%         - cell array B, the polynomial basis to be used 
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - vector tau, the vector of collocation points
%         - string cost_policy, the policy to be minimized
%         - string basis, the polynomial basis to be used
%         - string method, the parameter distribution to be used

% Outputs: - scalar r, the cost index to be optimized

function [r] = cost_function(mu, initial, final, measurements, n, x, B, cost_policy, basis, tau, sampling_distribution)
    % Minimize the control input
    P = reshape(x(1:end-2), [length(n), max(n)+1]);     % Control points
    tf = x(end-1);                                      % The final time of flight
    N = floor(x(end));                                  % The optimal number of revolutions

    % Boundary conditions
    P = boundary_conditions(tf, n, initial, final, N, P, B, basis);

    % Re-evaluate the measurement's epochs in case it is needed
    epochs = measurements(1,:)/tf; 
    switch (sampling_distribution)
        case 'Regularized'
            options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);
            [~, epochs] = ode45(@(t,s)Sundman_transformation(basis,n,P,t,s), epochs, 0, options);
            epochs = epochs.';
        case 'Chebyshev'
            epochs = 2*epochs-1;
        case 'Legendre'
            epochs = 2*epochs-1;
        case 'Laguerre'
            epochs = 2*epochs-1;
    end

    % Cost function
    switch (cost_policy)
        case 'Least Squares'
            % State evolution
            B = state_basis(n, epochs, basis);
            C = evaluate_state(P,B,n);
        
            % Compute the residuals 
            M = cylindrical2cartesian(C(1:3,:),true);
            M(1:3,:) = M(1:3,:)./sqrt(M(1,:).^2+M(2,:).^2+M(3,:).^2);
            e = measurements(2:4,:)-M(1:3,:);
            r = sum(dot(e,e,1));  

        case 'Dynamics residual'
            % State evolution
            C = evaluate_state(P,B,n);

            % Control input
            u = acceleration_control(mu,C,tf,sampling_distribution);        

            % Control cost
            switch (sampling_distribution)
                case 'Regularized'
                    r = sqrt(C(1,:).^2+C(3,:).^2);                   % Radial evolution
                    a = sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2);         % Non-dimensional acceleration
                    a = a./r;                                        % Dimensional acceleration
                otherwise
                    a = sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2);         % Dimensional acceleration
            end
            
            % Cost function
            r = trapz(tau,a); 

        otherwise 
            error('No valid cost function was selected');
    end
end

%% Auxiliary function
% Compute the Sundman transformation 
function [ds] = Sundman_transformation(basis, n, P, t, s)
    B = state_basis(n, s, basis);
    C = evaluate_state(P,B,n);
    ds = (C(1,:).^2+C(3,:).^2).^(-1/2);
end