%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    % Differential time law
    l = 1 + s(1,:) .* cos(t) + s(2,:) .* sin(t);
    gamma = s(3,:).^3 .* l.^2;

    % Mayer and Lagrange terms
    M = -0 * s(3,end)^2; 
    L = dot(u, u, 1) ./ gamma;
end