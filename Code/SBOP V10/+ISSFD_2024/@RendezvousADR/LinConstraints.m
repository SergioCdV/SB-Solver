%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - class Problem, defining the problem at hands

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq    

function [A, b, Aeq, beq] = LinConstraints(obj, params, beta, P)
    % Constants 
    Dim = 2+size(beta,1)+size(P,1)*size(P,2);    % Total dimension of the optimization variables

    % Linear inequalities
    A = zeros(Dim);
    A(1,end-size(beta,1)-1) = 1;        % The initial time must be smaller than the final time (the independent variable is monotone)
    A(1,end-size(beta,1)) = -1;
    b = zeros(Dim,1);
    
    % Linear constraints
    Aeq = zeros(Dim);
    Aeq(1,end-size(beta,1)-1) = 1;      % Initial true anomaly constraints
    beq = zeros(Dim,1);
    beq(1,1) = params(5);

    A(2,end-size(beta,1)) = 1;          % Final true anomaly inequality
    b(2,1) = params(6);
end