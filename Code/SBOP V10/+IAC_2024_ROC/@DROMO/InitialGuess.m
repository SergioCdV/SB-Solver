%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Initial guess function %% 
% Function implementation of the a warming up initial guess if available

function [beta, t0, tf] = InitialGuess(obj, params, initial, final)    
    % New initial TOF (anomaly)
    t0 = 0;
    tf = params(5) + 2 * pi * params(6) + atan2(final(2), final(1));
    beta = final(1:2);
end