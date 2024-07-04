%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    % Sundman transformation 
    mu = params(1);                                                 % Gravitational parameter

    l = t(1,:) + s(6,:);
    w = 1 + s(2,:) .* cos(l) + s(3,:) .* sin(l);

    gamma = sqrt(mu * s(1,:)) .* (w ./ s(1,:)).^2;

    % Minimum energy cost function 
    M = 0;%-s(1,end); 
    L = sqrt( dot(u,u,1) ) .* gamma;%zeros(1,size(u,2)); %0;
end