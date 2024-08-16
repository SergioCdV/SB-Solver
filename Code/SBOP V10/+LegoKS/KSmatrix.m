%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 17/11/22

%% KS matrix %%
% Function to compute the operator to transform from the position space to
% the u-space

% Inputs: - vector u, the state variable

% Outputs: - array L, the KS operator

function [L] = KSmatrix(u)
    % Compute the L operator 
    L = [u(1) -u(2) -u(3) u(4); u(2) u(1) -u(4) -u(3); u(3) u(4) u(1) u(2); u(4) -u(3) u(2) -u(1)];
end