%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Initial guess function %% 
% Function implementation of the a warming up initial guess if available

function [beta, t0, tf] = InitialGuess(obj, params, initial, final)    
    % New initial TOF
    t0 = params(end-2);
    tf = params(end-1) + 2 * pi * params(end);
    beta = [0; 0];
end