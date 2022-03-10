%% Project: 
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - vector n, the approximation degree to each position coordinate
%         - array P, the set of control points to estimate the position vector 
%         - array B, the polynomial basis in use in the approximation

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(tf, tau, initial, final, n, m, x, B)
    % Extract the optimization variables
    P = reshape(x(1:end-m-1), [1, max(n)+1]);
    u = reshape(x(end-m:end-1), [1 m]);
    C = evaluate_state(P,B,n);

    % Non-linear inequality
    c = [];

    % Boundary conditions
    ceq = [C(1,1)-initial.'; C(1,end)-final.'];

    C = C(:,1:end);
    u = u(1,1:end);

    % Dimensionalising 
    C(2,:) = C(2,:) / (tf * x(end));

    % Dynamic constraints  
    D = C(2,:)-u+C(1,:).^3;
    ceq = [ceq; reshape(D, [size(D,1)*size(D,2) 1])];
end