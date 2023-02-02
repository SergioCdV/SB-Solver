%% Project: 
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - vector n, the approximation degree to each position coordinate
%         - array P, the set of control points to estimate the position vector 
%         - array B, the polynomial basis in use in the approximation

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(tf, tau, initial, n, m, x, B)
    % Extract the optimization variables
    P = reshape(x(1:end-m), [1, max(n)+1]);
    u = reshape(x(end-m+1:end), [1 m]);
    C = evaluate_state(P,B,n);

    % Non-linear inequality
    c = [];

    % Boundary conditions
    ceq = C(1,1)-initial(1);

    % Dynamic constraints  
    C = C(:,2:end-1);
    u = u(2:end-1);
    D = C(2,:)-u;
    ceq = [ceq; reshape(D, [size(D,1)*size(D,2) 1])];
end