%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Initial guess function %% 
% Function implementation of the a warming up initial guess if available

function [beta, t0, tf] = InitialGuess(obj, params, initial, final)    
    % Initial guess
    t0 = params(1);            % Initial clock
    tf = params(2);            % Final clock

    if (length(params) > 58)
        beta = params(59:64).';    % Final joints states
    else
        beta = initial;            % Final joints states
    end
end