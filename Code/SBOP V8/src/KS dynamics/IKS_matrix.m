%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 17/11/22

%% Inverse KS matrix %%
% Function to compute the inverse of the operator to transform from the position space to
% the u-space

% Inputs: - vector u, the state variable

% Outputs: - array L, the inverse of the KS operator

function [L] = IKS_matrix(u)
    % Compute the L operator 
    L = L_operator(u); 
    
    % Compute the inverse of the L operator
    L = L.'/norm(u)^2;
end