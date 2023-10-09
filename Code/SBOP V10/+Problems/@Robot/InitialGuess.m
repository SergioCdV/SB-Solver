%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Initial guess function %% 
% Function implementation of the a warming up initial guess if available

function [beta, t0, tf] = InitialGuess(obj, params, initial, final)    
    % Initial guess
    t0 = params(26);                    % Initial true anomaly [rad]
    tf = params(27);                    % Final true anomaly [rad]
    beta = [100 100 100 0 0 0].';
end