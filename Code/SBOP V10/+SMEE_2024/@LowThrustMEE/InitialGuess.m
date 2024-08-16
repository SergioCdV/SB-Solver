%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Initial guess function %% 
% Function implementation of the a warming up initial guess if available

function [beta, t0, tf] = InitialGuess(obj, params, initial, final)    
    % New initial TOF (anomaly)
    t0 = 0;
    tf = params(6);
    beta = [initial(1); 0];
end