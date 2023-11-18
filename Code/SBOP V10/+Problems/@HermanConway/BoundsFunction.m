%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 10/11/2023

%% Uppe an lower bounds function %% 
% Function implementation the definition of the upper na dlowe bounds for
% the problem

function [LB, UB] = BoundsFunction(obj)
    % Upper and lower bounds for the problem first order state vector, initial time, final time and parameters
    LB = [-Inf * ones(1,2) 0 0  -Inf * ones(1,4)];
    UB = [ Inf * ones(1,2) 1 100 Inf * ones(1,4)];
end