%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    % Differential time law 
    omega = params(23)^2 / params(25)^3;       % True orbit mean motion
    k = 1 + params(24) * cos(t(1,:));          % Transformation parameter
    dnu = omega .* k;                          % Differential time law (Kepler's second law)

    % Mayer and Lagrange terms
    M = 0; 
    L = dot(u, u, 1) ./ dnu;
end