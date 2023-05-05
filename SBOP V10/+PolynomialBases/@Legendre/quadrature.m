%% Project: SBOPT %%
% Date: 05/05/23

%% Quadrature scheme %%
% Function to compute the quadrature scheme weights and associated
% integration domain

% Inputs: 

% Outputs: - vector tau, the domain of integration
%          - vector W, the quadrature weights 
%          - scalar J, the Jacobian of the domain transformation

function [tau, W, J, f, D] = quadrature(n, m)
    % Sanity check on the quadrature number of points 
    if ((m-1)/2 <= max(n))
        m = 2*max(n)+1;
    end

    % Jacobian domain transformation 
    J = 0.5;

    % Guass-Legendre Quadrature weights
    [tau, W, D] = weights(m);

    % Domain transformation 
    f = @(t0, tf, tau)[(tf - t0) * 0.5 * (1+tau); (tf-t0) * J * ones(1,length(tau))];       
end