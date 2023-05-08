%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, s)
    % Compute the control vector as a dynamics residual
    theta = atan2(s(4,:), s(3,:));
    u = s(3:4,:)-params(1) * [cos(theta); sin(theta)];
end