%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    % Sundman transformation
    dtheta = dot(s(1:4,:), s(1:4,:), 1);  

    % Cost function
    M = 0;                                            % Mayer term
    L = dot(u(1:4,:), u(1:4,:), 1) .* dtheta;         % Minimum time transfer
end