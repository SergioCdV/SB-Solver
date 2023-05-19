%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, s)
    % Constants 
    mu = params(1); 
    
    % Linear terms of the equations of motion
    f = -mu * s(5,:)/4 .* s(1:4,:);                 % Acceleration vector
    a = s(11:14,:);                                 % Inertial acceleration field

    % Compute the control vector as a dynamics residual
    u = 2 * (a-f);
end