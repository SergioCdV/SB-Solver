%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    % Differential time law
    l = 2*s(1,:)-s(1,:).^2-s(2,:).^2;
    gamma = (s(3,:).*l).^(3/2).*s(1,:);

    % Mayer and Lagrange terms
    M = 0; 
    L = dot(u,u,1).*gamma;
end