%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    % Minimum fuel 
    M = 0; 

    S = 1+s(2,:).*cos(t)+s(3,:).*sin(t);
    gamma = s(3,:).^3.*S.^2;
    L = dot(u,u,1)./gamma;
end