%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Boundary Conditions function %% 
% Function implementation of the boundary conditions definition

function [s0, sf] = BoundaryConditions(obj, initial, final, beta, t0, tf)
    s0 = initial; 
    sf = [beta(1:6); zeros(3,1); beta(7:9)];
end