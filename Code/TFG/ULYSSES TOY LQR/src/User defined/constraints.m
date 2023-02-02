%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - vector initial, the initial boundary conditions of the
%           trajectory 
%         - vector final, the initial boundary conditions of the
%           trajectory 
%         - matrix D, the differentiation matrix
%         - vector m, the number of nodes
%         - vector tau, the normalized independent variable
%         - vector x, the degree of freedom to be optimized

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(initial, final, tf, D, m, tau, x)
    % Extract the optimization variables
    C = reshape(x(1:2*(m+1)), [2 length(tau)]);         % State points
    u = reshape(x(end-m:end), [1 length(tau)]);         % Control points

    % Boundary equalities 
    bcs = [C(1:2,1)-initial; C(1:2,end)-final];

    % Path constraints 

    % Dynamics 
    res = acceleration_control(D, C, u, tf, tau);

    % Final constraints
    ceq = [bcs; reshape(res, [], 1)];              % Boundary and dynamics constraints
    c = [];                                        % Path inequalities
end