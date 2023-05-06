%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, s)
    % Constants 
    g = params(1); 
  
    % Compute the control vector as a dynamics residual
    u = s(3,:)-g*ones(1,length(t));
end