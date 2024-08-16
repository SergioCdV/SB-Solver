%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    % Sundman transformation 
    mu = params(1);                                         % Gravitational parameter
    
    l = t(1,:);
    w = 1 + s(2,:) .* cos(l) + s(3,:) .* sin(l);
    k = s(4,:) .* sin(l) - s(5,:) .* cos(l);
    gamma = sqrt(mu * s(1,:)) .* (w ./ s(1,:)).^2 + sqrt(s(1,:) / mu) .* k ./ w .* u(3,:);

    % Minimum energy cost function 
%     a = s(1,end) / (1 - s(2,end)^2 - s(3,end)^2);
    M = 0;%mu / (2 * a); 
    L = sqrt( dot(u,u,1) ) .* gamma; 
end 