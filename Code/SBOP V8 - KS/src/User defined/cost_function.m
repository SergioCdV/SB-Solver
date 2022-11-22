%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 31/01/22

%% Cost function %%
% Function to compute the cost function to be minimized

% Inputs: - string cost, the cost function to be minimized
%         - scalar mu, the gravitational parameter of the central body
%         - vector initial, the initial boundary conditions of the
%           trajectory
%         - vector final, the initial boundary conditions of the
%           trajectory
%         - cell array B, the polynomial basis to be used
%         - string basis, the polynomial basis to be used
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - vector tau, the vector of collocation points
%         - vector W, the quadrature weights
%         - vector x, the degree of freedom to be optimized
%         - string dynamics, the independent variable parametrization to be
%           used

% Outputs: - scalar r, the cost index to be optimized

function [r] = cost_function(cost, mu, initial, final, B, basis, n, tau, W, x)
    % Optimization variables
    tf = x(end-1);              % Final time of flight

    % Compute the state variable
    P = reshape(x(1:end-2), [length(n), max(n)+1]);                             % Control points
    P = boundary_conditions(tf, n, initial, final, P, B, basis);                % Boundary conditions control points
    C = evaluate_state(P,B,n);                                                  % State evolution
    
    % Compute the cost function
    switch (cost)
        case 'Minimum time'
            a = dot(C(1:4,:), C(1:4,:), 1);                                     % Time transformation
    
        case 'Minimum fuel'
            [u, ~, ~] = acceleration_control(mu, C, tf);                        % Control vector
            u = u / tf^2;                                                       % Normalized control vector
        
            a = sqrt(dot(u,u,1)).*dot(C(1:4,:), C(1:4,:), 1);                                               % Non-dimensional acceleration norm
    
        otherwise
            error('No valid cost function was selected to be minimized');
    end

    % Cost function
    if (isempty(W))
        r = tf*trapz(tau,a);
    elseif (length(W) ~= length(tau))
        r = 0; 
        for i = 1:floor(length(tau)/length(W))
            r = r + tf*dot(W,a(1+length(W)*(i-1):length(W)*i));
        end
    else
        r = tf*dot(W,a);
    end
end