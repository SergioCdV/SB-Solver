%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, s)
    % Compute the control vector as a dynamics residual
    u = (s(2,:) - 2 * s(1,:))./(2 * sqrt(s(1,:)));
end