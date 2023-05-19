%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Boundary Conditions function %% 
% Function implementation of the boundary conditions definition

function [s0, sf] = BoundaryConditions(obj, initial, final, beta, t0, tf)
    s0 = initial; 
    sf = final;
    thetaf = beta(1);
    R = [cos(thetaf) 0 0 -sin(thetaf); 0 cos(thetaf) sin(thetaf) 0; 0 -sin(thetaf) cos(thetaf) 0; sin(thetaf) 0 0 cos(thetaf)];
    sf([1:4 6:9]) = blkdiag(R,R) * final([1:4 6:9]);
end