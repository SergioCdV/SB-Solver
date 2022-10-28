%% Project: Shape-based attitude planning %%
% Date: 01/02/22

%% Hat map %%
% This file contains the function to implement the hat map of a R3 vector
% to the Lie algebra of SO3

% Inputs: - vector v, the vector to be mapped.

% Ouputs: - matrix S, the hat map of v.

% New version updates: 

function [S] = hat_map(v)
    %Hat map of the vector v
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end