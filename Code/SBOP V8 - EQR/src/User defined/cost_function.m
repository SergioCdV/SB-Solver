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
%         - vector uprev, the previously converged control profile
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

function [r] = cost_function(cost, mu, initial, final, uprev, B, basis, n, L, tau, W, x)
    % State evolution
    P = reshape(x(1:end-2), [length(n), max(n)+1]);                             % Control points
    thetaf = x(end-1);                                                          % Final insertion phase
    P = boundary_conditions(n, initial(1:5), final(1:5), P, B, basis);          % Boundary conditions control points
    C = evaluate_state(P,B,n);                                                  % State evolution

    % Compute the longitude evolution 
    theta0 = initial(end)-thetaf*L(1);
    L = theta0+thetaf*L;

    % Control vector 
    [u, ~, ~] = acceleration_control(mu, C, L, uprev);                          % Control vector
    u = u/thetaf;

    switch (cost)
        case 'Minimum time'
            % Kinematic constraint 
            w = 1+C(2,:).*cos(L)+C(3,:).*sin(L);
            a = sqrt(mu*C(1,:)).*(w./C(1,:)).^2;
            for i = 1:size(C,2)
                B = control_input(mu, C(:,i)); 
                a(i) = 1/(a(i)+B(6,3)*u(3,i));
            end
    
        case 'Minimum fuel'
            a = sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2);                           % Non-dimensional acceleration norm
        
        otherwise
            error('No valid cost function was selected to be minimized');
    end

    % Cost function
    if (isempty(W))
        r = trapz(tau,a);
    elseif (length(W) ~= length(tau))
        r = 0; 
        for i = 1:floor(length(tau)/length(W))
            r = r + dot(W,a(1+length(W)*(i-1):length(W)*i));
        end
    else
        r = dot(W,a);
    end

end
