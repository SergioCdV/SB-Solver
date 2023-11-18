%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    M = params(1)/s(1,end) - 0.5 * (s(3,end)^2 + (s(1,end) * s(4,end))^2); 
    L = zeros(1,size(u,2));
end