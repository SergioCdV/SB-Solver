%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    % Parameters 
    mu = params(1); 

    % Maximum energy 
    a = s(1,end) / (1 - s(2,end)^2 - s(3,end)^2);
    M = mu / (2 * a);
    L = zeros(1,size(u,2));
end