%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    % Sundman transformation 
%     mu = params(1);     % Gravitational parameter
%     w = 1 + s(2,:) .* cos(t(1,:)) + s(3,:) .* sin(t(1,:));
%     k = s(4,:) .* sin(t(1,:)) - s(5,:) .* cos(t(1,:));
% 
%     gamma = sqrt(mu * s(1,:)) .* (w ./ s(1,:)).^2 + sqrt(s(1,:) / mu) .* k ./ w .* u(3,:);

    % Minimum energy cost function 
    M = 0; 
    L = sqrt( dot(u,u,1) );
end