%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, tau, s)
    % Compute the control vector as a dynamics residual
    A = reshape(params(2:end), 3, []);

    % Control law
    u = s(7:9,:) - A * s(1:6,:);
end