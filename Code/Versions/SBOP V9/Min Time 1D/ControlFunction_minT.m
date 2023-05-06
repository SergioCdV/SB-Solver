%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(params, beta, t0, tf, tau, s)
    % Compute the control vector as a dynamics residual
    u = s(3,:)+sin(s(3,:));
end