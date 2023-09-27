%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Boundary Conditions function %% 
% Function implementation of the boundary conditions definition

function [s0, sf] = BoundaryConditions(obj, initial, final, beta, t0, tf)
    s0 = initial; 
    sf = final;

    s0(5:8) = beta(2:5);
    sf(1:8) = beta(6:13);
end