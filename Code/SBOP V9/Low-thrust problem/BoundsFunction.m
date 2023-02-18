%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Uppe an lower bounds function %% 
% Function implementation the definition of the upper na dlowe bounds for
% the problem

function [LB, UB] = BoundsFunction()
    % Upper and lower bounds for the problem first order state vector, initial time, final time and parameters
    LB = [-10 -10 -10 0 0 0];
    UB = [10 10 10 1 1e2 2*pi*1e2];
end