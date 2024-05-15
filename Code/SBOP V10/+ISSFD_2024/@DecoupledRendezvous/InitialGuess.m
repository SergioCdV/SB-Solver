%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Initial guess function %% 
% Function implementation of the a warming up initial guess if available

function [beta, t0, tf] = InitialGuess(obj, params, initial, final)    
    % Initial guess
    t0 = params(5);                    % Initial true anomaly [rad]
    tf = params(6);                    % Final true anomaly [rad]

    if ( length(params) <= 26 )
        beta = zeros(3,1); 
    else
        beta = params(27:29)';
    end
end