%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    % Mayer and Lagrange terms
    M = 0 * s(1,end); 
    L = 0 * dot(s(7:12,:), s(7:12,:), 1);
end