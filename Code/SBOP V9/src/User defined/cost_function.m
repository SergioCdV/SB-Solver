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

function [r] = cost_function(Cost, initial, final, B, basis, n, tau, W, x)
    % Minimize a given control input
    P = reshape(x(1:end-3), [length(n), max(n)+1]);                             % Control points
    thetaf = x(end-1);                                                          % Final insertion phase
    P = boundary_conditions(tf, n, initial, final, thetaf, P, B, basis);        % Boundary conditions control points
    C = evaluate_state(P,B,n);                                                  % State evolution
        
    % Evaluate the cost function
    [L, M] = Cost.evaluate_cost(beta, t0, tf, C, u); 
    r = M; 

    if (isempty(W))
        r = r + (tf-t0) * trapz(tau,L);
    else
        r = r + (tf-t0) * dot(W,L);
    end
end
