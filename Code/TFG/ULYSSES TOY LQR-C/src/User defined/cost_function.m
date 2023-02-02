%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 31/01/22

%% Cost function %%
% Function to compute the cost function to be minimized

% Inputs: - scalar tf, the final time of flight
%         - vector m, the number of nodes
%         - vector tau, the vector of collocation points
%         - vector W, the quadrature weights
%         - vector x, the degree of freedom to be optimized

% Outputs: - scalar r, the cost index to be optimized

function [r] = cost_function(tf, m, tau, W, x)    
    % Minimize the control input
    u = reshape(x(end-m:end), [1 length(tau)]);     % Control points        
    a = 0.5*u(1,:).^2;                              % Non-dimensional acceleration norm

    % Cost function
    if (isempty(W))
        r = tf*trapz(tau,a);
    elseif (length(W) ~= length(tau))
        r = 0; 
        for i = 1:floor(length(tau)/length(W))
            r = r + tf*dot(W,a(1+length(W)*(i-1):length(W)*i));
        end
    else
        r = tf/2*dot(W,a);
    end
end
