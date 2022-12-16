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
    % Extract the optimization variables
    P = reshape(x(1:end-3), [length(n), max(n)+1]);     % Control points
    thetaf = x(end-2);                                  % Final fiber parameter
    sf = x(end-1);                                      % Final time of flight 

    R = [cos(thetaf) 0 0 -sin(thetaf); 0 cos(thetaf) sin(thetaf) 0; 0 -sin(thetaf) cos(thetaf) 0; sin(thetaf) 0 0 cos(thetaf)];
    final = final*blkdiag(R,R).';

    P = boundary_conditions(sf, n, initial, final, P, B, basis);                % Boundary conditions control points
    C = evaluate_state(P,B,n);                                                  % State evolution

    % Sundman transformation 
    r = dot(C(1:4,:),C(1:4,:),1);
    
    % Compute the cost function
    switch (cost)
        case 'Minimum time'
            a = r;                                                              % Time transformation
    
        case 'Minimum fuel'
            [u, ~, ~] = acceleration_control(mu, C, sf);                        % Control vector
            u = u(1:3,:) / sf^2;                                                % Normalized control vector
        
            a = sqrt(dot(u,u,1));                                               % Non-dimensional acceleration norm
    
        otherwise
            error('No valid cost function was selected to be minimized');
    end

    % Cost function
    if (isempty(W))
        r = sf*trapz(tau,a);
    elseif (length(W) ~= length(tau))
        r = 0; 
        for i = 1:floor(length(tau)/length(W))
            r = r + sf*dot(W,a(1+length(W)*(i-1):length(W)*i));
        end
    else
        r = sf*dot(W,a);
    end
end
