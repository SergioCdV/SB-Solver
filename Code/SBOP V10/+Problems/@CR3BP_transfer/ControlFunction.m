%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, tau, s)
    % Constants 
    mu = params(1);         % Gravitational parameter of the system 
    
    % Compute the Coriolis terms 
    Gamma = [-2 * s(5,:) - s(1,:); + 2 * s(4,:) - s(2,:); zeros(1,size(s,2))];

    % Compute the dipoles
    f = - (1-mu) * s(1:3,:) ./ sqrt( dot(s(1:3,:), s(1:3,:), 1) ).^3 - mu * s(1:3,:) ./ sqrt( dot(s(1:3,:), s(1:3,:), 1) ).^3;

    % Compute the control vector as a dynamics residual
    u(1:3,:) = s(7:9,:) + Gamma - f;
end