%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - scalar tf, the final time of flight
%         - vector initial, the initial boundary conditions of the
%           trajectory 
%         - vector final, the initial boundary conditions of the
%           trajectory 
%         - cell array B, the polynomial basis to be used
%         - string basis, the polynomial basis to be used
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - vector tau, the normalized independent variable
%         - vector x, the degree of freedom to be optimized

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(tf, initial, final, B, basis, n, tau, x)
    % Extract the optimization variables
    P = reshape(x, [length(n), max(n)+1]);                                      % Control points

    % Boundary conditions points
    P = boundary_conditions(tf, n, initial, final, P, B, basis);

    % Trajectory evolution
    C = evaluate_state(P,B,n);

    % Equalities constraints
    ceq = [];

    % Inequalities constraints
    c = [];
end