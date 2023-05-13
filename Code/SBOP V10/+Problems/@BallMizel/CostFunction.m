%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, tau, s, u)
    M = 0; 
    L = (tau.^2-s(1,:).^3).^2 .* u.^14 + params(1) * u.^2;
end