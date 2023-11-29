%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, tau, s)
    % Compute the control vector as a dynamics residual (linear acceleration, TH relative motion model)
    u(1:3,:) = s(7:9,:) - [2 * s(6,:); -s(2,:); 3*s(3,:)-2*s(4,:)];
end