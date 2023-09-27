%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    % Orbital dynamics
    COE = params(2:7);                                      % Orbit elements 
    v = (tf-t0) * 0.5 * (1+t) + COE(end);                   % Anomaly space 
    h = sqrt(params(1) * COE(1) * (1-COE(2)^2));            % Angular momentum 
    r = COE(1) * (1-COE(2)^2) ./ (1 + COE(2)*cos(v));       % Spacecraft position vector
    dtheta = h./r.^2;                                       % Angular velocity of the anomaly

    % Cost functions
    M = 0; 
    L = (0 * ones(1,size(u,2)) + 1 * dot(u,u,1)) ./ dtheta;
end