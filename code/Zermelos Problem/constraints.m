%% Project: 
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - vector n, the approximation degree to each position coordinate
%         - array P, the set of control points to estimate the position vector 
%         - array B, the polynomial basis in use in the approximation

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(V, w, initial, final, n, m, x, B)
    % Extract the optimization variables
    P = reshape(x(1:end-m), [2, max(n)+1]);

    % Non-linear inequality
    c = [];

    % Boundary conditions
    C = evaluate_state(P,B,n);
    ceq = [C(1:2,1)-initial(1:2); C(1:2,end)-final(1:2)];

    % Dynamic constraints 
    theta = x(end-m+2:end-1).';
    u = [cos(theta); sin(theta)];
    D = C([3 4], 2:end-1)-w-V*u; 
    ceq = [ceq; reshape(D, [size(D,1)*size(D,2) 1])];
end