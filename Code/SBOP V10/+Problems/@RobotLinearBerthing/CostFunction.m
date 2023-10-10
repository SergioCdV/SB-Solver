%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    % Differential time law 
    rho = 1 + params(5) * cos(t(1,:));                    % Transformation parameter
    Omega = params(4)^2/params(6)^3 .* rho.^2;            % True anomaly angular velocity [rad/s]

    % Mayer and Lagrange terms
    M = 0; 
    L = dot(u, u, 1) .* Omega;
end