%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    % Sundman transformation
    r = dot(s(1:4,:), s(1:4,:), 1);
    dtheta = r; 

    % Correct dimensions
    u = u ./ r.^2;

    % Cost function
    M = 0;% * -r(end);                      % Mayer term
    l = dot(u(1:4,:), u(1:4,:), 1);         % Minimum time transfer
    L = l .* dtheta;
end