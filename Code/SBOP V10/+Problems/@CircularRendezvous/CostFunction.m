%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    % Mayer and Lagrange terms
    M = 0;
    if (length(params) < 12)
        L = dot(u, u, 1);
    else
        L = dot(u, u, 1) + ones(1, size(u,2));
    end
end