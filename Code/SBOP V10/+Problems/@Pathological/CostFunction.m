%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    M = 0; 
    L = (u.^2 - s(1,:).^2) ;
end