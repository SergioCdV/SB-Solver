%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Boundary Conditions function %% 
% Function implementation of the boundary conditions definition

function [s0, sf] = BoundaryConditions(obj, initial, final, params, beta, t0, tf)
    % Initial fiber
    R = LegoKS.GroupAction( beta(end-1) );
    R = blkdiag(R,R);
    s0 = R * initial; 
    
    % Final fiber
%     final = beta(end-9:end-2);
    R = LegoKS.GroupAction( beta(end) );
    R = blkdiag(R,R);
    sf = R * final;
end