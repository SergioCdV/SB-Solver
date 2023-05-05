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
    if (m <= max(n))
        m = max(n)+1;
    end

    % Jacobian domain transformation 
    J = 1; 

    % Scaled Gauss-Legendre Quadrature weights
    [tau, W, D] = weights(m);
    W = W/2;
    tau = 0.5*(tau+1);

    % Domain transformation 
    f = @(t0, tf, tau)[(tf-t0) * tau; (tf-t0) * ones(1,length(tau))];        
end