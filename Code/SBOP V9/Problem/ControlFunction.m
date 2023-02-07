%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(m, beta, t0, tf, tau, s)
    u = zeros(m,size(s,2));
end